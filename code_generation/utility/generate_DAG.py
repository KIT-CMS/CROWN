import re
import json
import os
from code_generation.quantity import QuantityGroup


class GraphParser:
    def __init__(self, config, nanoAOD_inputs, DAG_dir=None, debugging=False):
        self.config = config
        self.nanoAOD_inputs = nanoAOD_inputs
        self.debugging = debugging
        self.nodes = {}
        self.edges = {}
        self.inputs = {}
        self.outputs = {}
        self.vec_output_mappings = {}
        self.in_file_name = "IN_File"
        self.out_file_name = "OUT_File"
        self.DAG_data = {}
        self.DAG_dir = DAG_dir

    def save_graph(self, scope, name):
        path = f"{name}_{self.config.era}_{self.config.sample}_{scope}.json"
        if self.DAG_dir:
            path = os.path.join(self.DAG_dir, path)
        scope_DAG_data = self.DAG_data[scope]
        if scope != "global":
            scope_DAG_data += self.DAG_data["global"]
        with open(path, "w") as f:
            json.dump(scope_DAG_data, f, indent=4, sort_keys=True)
        self.update_DAG_file_list(os.path.join(self.DAG_dir, "DAG_files.json"),self.config.era, self.config.sample, scope)


    def update_DAG_file_list(self, config_path, new_era, new_sample, new_scope):
        """
        Updates the config JSON with new values, ensuring no duplicates.
        """
        # Initialize with empty structure if file doesn't exist
        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                try:
                    data = json.load(f)
                except json.JSONDecodeError:
                    data = {"era": [], "sample": [], "scope": []}
        else:
            data = {"era": [], "sample": [], "scope": []}

        # Add new values and ensure uniqueness using sets
        # Then convert back to sorted lists for a clean JSON output
        data["era"] = sorted(list(set(data["era"] + [str(new_era)])))
        data["sample"] = sorted(list(set(data["sample"] + [str(new_sample)])))
        data["scope"] = sorted(list(set(data["scope"] + [str(new_scope)])))

        # Write back to the file
        with open(config_path, 'w') as f:
            json.dump(data, f, indent=4)

        print(f"Updated {config_path} with {new_era}/{new_sample}/{new_scope}")

    def add_output(self, output_name, output_node, scope):
        if not self.outputs.get(scope):
            self.outputs[scope] = {}
        if not self.outputs[scope].get(output_name):
            self.outputs[scope][output_name] = []
        self.outputs[scope][output_name] += [f"{scope}_{output_node}"]

    def add_input(self, input_name, input_node, scope):
        scoped_input = f"{scope}_{input_node}"
        if not self.inputs.get(scope):
            self.inputs[scope] = {}
        if not self.inputs[scope].get(scoped_input):
            self.inputs[scope][scoped_input] = []
        self.inputs[scope][scoped_input] += [input_name]

    def add_node(self, id_name, name=None, scope="global", parent=None, node_type=None):
        if not name:
            name = id_name
        add_data = {"data": {"id": f"{scope}_{id_name}", "label": name}}
        if parent:
            add_data["data"]["parent"] = f"{scope}_{parent}"
        if node_type:
            add_data["data"]["type"] = node_type
        if not self.nodes.get(scope):
            self.nodes[scope] = []
        self.nodes[scope] += [add_data]

    def add_edge(self, source, target, scope, name):
        if not self.edges.get(scope):
            self.edges[scope] = []
        self.edges[scope] += [{"data": {"source": source, "target": target, "label": name}}]

    def parse_from_call(self, producer, parent, scope, align=""):
        group_id = f"{producer.name}_v"
        call = producer.call
        pattern = r'event::filter::Flag\({df}, "{(.*)}", "{(.*)}"\)'
        match = re.search(pattern, call)
        if match:
            input_vec_config = match.group(1)
        else:
            raise Exception(f"Input vector config could not be parsed from {call}")
        if not input_vec_config in producer.vec_configs:
            raise Exception(
                f"Input name from {call} not in producer vector configs {producer.vec_configs}"
            )
        else:
            vec_input_index = producer.vec_configs.index(input_vec_config)

        self.add_node(
            id_name=group_id,
            name=producer.name,
            scope=scope,
            parent=parent,
            node_type="vector",
        )
        vec_configs = [
            self.config.config_parameters[scope][c] for c in producer.vec_configs
        ]

        for i_c, c in enumerate(zip(*vec_configs)):
            input_name = c[vec_input_index]
            vector_id = f"{group_id}_{i_c}"
            vector_name = f"{producer.name}_{i_c}"
            if self.debugging:
                print(align + "    " + f"Adding Producer: {vector_name}")
            self.add_node(
                id_name=vector_id, name=vector_name, scope=scope, parent=group_id
            )
            self.add_input(input_name, vector_id, scope)
            if isinstance(producer.output, list):
                raise NotImplementedError

    def generate_graph(self):

        self.add_node(id_name=self.in_file_name, node_type="I/O")

        for scope in self.config.scopes:
            if self.debugging:
                print()
                print(f"For scope {scope}:")
            self.add_node(
                id_name="scope", name=f"{scope} scope", scope=scope, node_type="scope"
            )
            for p in self.config.producers[scope]:
                self.parse_producer(
                    producer=p, parent="scope", scope=scope, align="    "
                )
            for output in self.config.outputs[scope]:
                if isinstance(output, QuantityGroup):
                    vec_config = self.config.config_parameters[scope][output.vec_config]
                    vec_output_name = self.vec_output_mappings[output.name]
                    for o in vec_config:
                        self.add_input(
                            o[vec_output_name], f"{self.out_file_name}", scope
                        )
                else:
                    self.add_input(output.name, f"{self.out_file_name}", scope)
            self.add_node(
                id_name=self.out_file_name,
                name=f"{scope}_{self.out_file_name}",
                scope=scope,
                node_type="I/O",
            )

        # check outputs for ambigous assignment
        for scope in self.config.scopes:
            all_keys = list(self.outputs["global"].keys()) + list(
                self.outputs[scope].keys()
            )
            for key in all_keys:
                total_output_nodes = self.outputs["global"].get(key, []) + self.outputs[
                    scope
                ].get(key, [])
                if len(set(total_output_nodes)) != 1:
                    raise Exception(
                        f"Output {output} has multiple origins: {total_output_nodes}"
                    )
        self.assemble_edges()
        for scope in self.config.scopes:
            self.DAG_data[scope] = self.nodes[scope] + self.edges[scope]

    def assemble_edges(self):
        for scope in self.config.scopes + ["global"]:
            for k, i in self.inputs[scope].items():
                compose = {}
                for i2 in i:
                    if self.outputs["global"].get(i2):
                        source = self.outputs["global"][i2][0]
                    elif self.outputs[scope].get(i2):
                        source = self.outputs[scope][i2][0]
                    elif i2 in self.nanoAOD_inputs:
                        source = f"global_{self.in_file_name}"
                    else:
                        raise Exception("WTF is wrong with this?")
                    if not compose.get(source):
                        compose[source] = []
                    compose[source] += [i2]
                for o, n in compose.items():
                    composite_name = ", ".join(sorted(n))
                    if self.debugging:
                        print(f"Adding Edge {composite_name} from {o} to {k}")
                    self.add_edge(source=o, target=k, scope=scope, name=composite_name)

    def parse_VectorProducer(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding VectorProducer: {producer.name}")
        print(
            f"!!!{producer} is a lagacy producer and should be replaced with ExtendedVectorProducer!!!"
        )
        if producer.call.startswith("event::filter::Flag"):
            self.parse_from_call(
                producer=producer, parent=parent, scope=scope, align=align
            )
        else:
            raise NotImplementedError(
                f"The call {producer.call} does not have legacy support."
            )

    def parse_ExtendedVectorProducer(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding VectorProducer: {producer.name}")
        group_id = f"{producer.name}_v"
        self.add_node(
            id_name=group_id,
            name=producer.name,
            scope=scope,
            parent=parent,
            node_type="vector",
        )
        vec_config = self.config.config_parameters[scope][producer.vec_config]
        for i_c, c in enumerate(vec_config):
            vector_id = f"{group_id}_{i_c}"
            vector_name = f"{producer.name}_{i_c}"
            if self.debugging:
                print(align + "    " + f"Adding Producer: {vector_name}")
            self.add_node(
                id_name=vector_id, name=vector_name, scope=scope, parent=group_id
            )
            if isinstance(producer.input[scope], list):
                for n in producer.input[scope]:
                    self.add_input(n.name, vector_id, scope)
                    if self.debugging:
                        print(align + "        " + f"Adding Input: {n.name}")
            if isinstance(producer.output, list):
                vec_output = c[producer.outputname]
                self.vec_output_mappings[producer.output_group.name] = (
                    producer.outputname
                )
                self.add_output(vec_output, vector_id, scope)
                if self.debugging:
                    print(align + "        " + f"Adding Output: {vec_output}")

    def parse_ProducerGroup(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding ProducerGroup: {producer.name}")
        group_id = f"{producer.name}_g"
        self.add_node(id_name=group_id, name=producer.name, scope=scope, parent=parent)
        for p in producer.producers[scope]:
            self.parse_producer(
                producer=p, parent=group_id, scope=scope, align=align + "    "
            )

    def parse_Filter(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding Filter: {producer.name}")
        group_id = f"{producer.name}_g"
        self.add_node(id_name=group_id, name=producer.name, scope=scope, parent=parent)
        for p in producer.producers[scope]:
            self.parse_producer(
                producer=p, parent=group_id, scope=scope, align=align + "    "
            )

    def parse_BaseFilter(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding Filter: {producer.name}")
        self.add_node(id_name=producer.name, scope=scope, parent=parent)
        for n in producer.input[scope]:
            self.add_input(n.name, producer.name, scope)
            if self.debugging:
                print(align + "    " + f"Adding Input: {n.name}")

    def parse_Producer(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding Producer: {producer.name}")
        self.add_node(id_name=producer.name, scope=scope, parent=parent)
        if isinstance(producer.input[scope], list):
            for n in producer.input[scope]:
                self.add_input(n.name, producer.name, scope)
                if self.debugging:
                    print(align + "    " + f"Adding Input: {n.name}")
        if isinstance(producer.output, list):
            for n in producer.output:
                self.add_output(n.name, producer.name, scope)
                if self.debugging:
                    print(align + "    " + f"Adding Output: {n.name}")

    def parse_producer(self, producer, parent, scope, align=""):
        if producer.__class__.__name__ == "VectorProducer":
            self.parse_VectorProducer(producer, parent, scope, align)
        elif producer.__class__.__name__ == "ExtendedVectorProducer":
            self.parse_ExtendedVectorProducer(producer, parent, scope, align)
        elif producer.__class__.__name__ == "ProducerGroup":
            self.parse_ProducerGroup(producer, parent, scope, align)
        elif producer.__class__.__name__ == "Filter":
            self.parse_Filter(producer, parent, scope, align)
        elif producer.__class__.__name__ == "BaseFilter":
            self.parse_BaseFilter(producer, parent, scope, align)
        elif producer.__class__.__name__ == "Producer":
            self.parse_Producer(producer, parent, scope, align)
        else:
            print(producer)
            raise Exception(f"Unknown Producer class {producer.__class__.__name__}")


def create_graph(configuration, nanoAOD_inputs, DAG_dir, json_name):
    graph = GraphParser(configuration, nanoAOD_inputs, DAG_dir, debugging=False)
    graph.generate_graph()
    for scope in graph.config.scopes:
        graph.save_graph(scope, json_name)
