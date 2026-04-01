import re
import json
import os
from collections import defaultdict
from code_generation.quantity import QuantityGroup

class GraphParser:
    """
    Parses configuration data to generate Directed Acyclic Graphs (DAGs)
    representing the dependencies and flow of producers, filters, and I/O.
    """
    def __init__(self, config, nanoAOD_inputs, DAG_dir=None, debugging=False):
        self.config = config
        self.nanoAOD_inputs = nanoAOD_inputs
        self.debugging = debugging
        self.nodes = defaultdict(list)
        self.edges = defaultdict(list)
        self.inputs = defaultdict(lambda: defaultdict(list))
        self.outputs = defaultdict(lambda: defaultdict(list))
        self.vec_output_mappings = {}
        self.in_file_name = "IN_File"
        self.out_file_name = "OUT_File"
        self.DAG_data = {}
        self.DAG_dir = DAG_dir

    def save_graph(self, scope, name):
        """
        Compiles the global and scoped DAG data and saves it to a JSON file.
        Also triggers an update to the master DAG file list.
        """
        path = f"{name}_{self.config.era}_{self.config.sample}_{scope}.json"
        if self.DAG_dir:
            path = os.path.join(self.DAG_dir, path)

        scope_DAG_data = self.DAG_data[scope]
        if scope != "global":
            scope_DAG_data += self.DAG_data["global"]

        with open(path, "w") as f:
            json.dump(scope_DAG_data, f, indent=4, sort_keys=True)

        self.update_DAG_file_list(
            os.path.join(self.DAG_dir, "DAG_files.json"),
            self.config.era,
            self.config.sample,
            scope
        )

    def update_DAG_file_list(self, config_path, new_era, new_sample, new_scope):
        """
        Updates the master JSON tracking file to ensure new eras, samples,
        and scopes are registered without duplicates.
        """
        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                try:
                    data = json.load(f)
                except json.JSONDecodeError:
                    data = {"era": [], "sample": [], "scope": []}
        else:
            data = {"era": [], "sample": [], "scope": []}

        data["era"] = sorted(list(set(data["era"] + [str(new_era)])))
        data["sample"] = sorted(list(set(data["sample"] + [str(new_sample)])))
        data["scope"] = sorted(list(set(data["scope"] + [str(new_scope)])))

        with open(config_path, 'w') as f:
            json.dump(data, f, indent=4)

        print(f"Updated {config_path} with {new_era}/{new_sample}/{new_scope}")

    def add_output(self, output_name, output_node, scope):
        """Registers a node as the creator of a specific output."""
        self.outputs[scope][output_name].append(f"{scope}_{output_node}")

    def add_input(self, input_name, input_node, scope):
        """Registers a node's requirement for a specific input."""
        scoped_input = f"{scope}_{input_node}"
        self.inputs[scope][scoped_input].append(input_name)

    def add_node(self, id_name, name=None, scope="global", parent=None, node_type=None):
        """Creates a standardized node object and appends it to the graph's node list."""
        if not name:
            name = id_name

        add_data = {"data": {"id": f"{scope}_{id_name}", "label": name}}
        if parent:
            add_data["data"]["parent"] = f"{scope}_{parent}"
        if node_type:
            add_data["data"]["type"] = node_type

        self.nodes[scope].append(add_data)

    def add_edge(self, source, target, scope, name):
        """Creates a directional edge between a source node and a target node."""
        self.edges[scope].append({"data": {"source": source, "target": target, "label": name}})

    def parse_from_call(self, producer, parent, scope, align=""):
        """
        Uses regex to extract vector configurations from a legacy producer's string call
        and maps them into graph nodes and inputs.
        """
        group_id = f"{producer.name}_v"
        call = producer.call
        pattern = r'event::filter::Flag\({df}, "{(.*)}", "{(.*)}"\)'
        match = re.search(pattern, call)

        if match:
            input_vec_config = match.group(1)
        else:
            raise ValueError(f"Input vector config could not be parsed from {call}")

        if input_vec_config not in producer.vec_configs:
            raise ValueError(f"Input name from {call} not in producer vector configs {producer.vec_configs}")

        vec_input_index = producer.vec_configs.index(input_vec_config)

        self.add_node(id_name=group_id, name=producer.name, scope=scope, parent=parent, node_type="vector")
        vec_configs = [self.config.config_parameters[scope][c] for c in producer.vec_configs]

        for i_c, c in enumerate(zip(*vec_configs)):
            input_name = c[vec_input_index]
            vector_id = f"{group_id}_{i_c}"
            vector_name = f"{producer.name}_{i_c}"

            if self.debugging:
                print(align + "    " + f"Adding Producer: {vector_name}")

            self.add_node(id_name=vector_id, name=vector_name, scope=scope, parent=group_id)
            self.add_input(input_name, vector_id, scope)

            if isinstance(producer.output, list):
                raise NotImplementedError("List outputs for legacy parsed calls are not supported.")

    def generate_graph(self):
        """
        The main orchestrator. Iterates through all scopes of the configuration,
        generates the inputs/outputs/nodes, verifies output ambiguity, and creates edges.
        """
        self.add_node(id_name=self.in_file_name, node_type="I/O")

        for scope in self.config.scopes:
            if self.debugging:
                print(f"\nFor scope {scope}:")
            self.add_node(id_name="scope", name=f"{scope} scope", scope=scope, node_type="scope")

            # Parse all producers in this scope
            for p in self.config.producers[scope]:
                self.parse_producer(producer=p, parent="scope", scope=scope, align="    ")

            # Parse outputs and connect to out_file_name
            for output in self.config.outputs[scope]:
                if isinstance(output, QuantityGroup):
                    vec_config = self.config.config_parameters[scope][output.vec_config]
                    vec_output_name = self.vec_output_mappings[output.name]
                    for o in vec_config:
                        self.add_input(o[vec_output_name], f"{self.out_file_name}", scope)
                else:
                    self.add_input(output.name, f"{self.out_file_name}", scope)

            self.add_node(
                id_name=self.out_file_name,
                name=f"{scope}_{self.out_file_name}",
                scope=scope,
                node_type="I/O",
            )

        # Check outputs for ambiguous assignments (multiple origins for the same output)
        for scope in self.config.scopes:
            all_keys = list(self.outputs["global"].keys()) + list(self.outputs[scope].keys())
            for key in all_keys:
                total_output_nodes = self.outputs["global"].get(key, []) + self.outputs[scope].get(key, [])
                if len(set(total_output_nodes)) != 1:
                    raise ValueError(f"Output {key} has multiple origins: {total_output_nodes}")

        self.assemble_edges()

        # Compile final DAG data
        for scope in self.config.scopes:
            self.DAG_data[scope] = self.nodes[scope] + self.edges[scope]

    def assemble_edges(self):
        """
        Resolves the inputs vs. outputs mappings to draw the actual connecting
        lines (edges) between nodes in the DAG.
        """
        for scope in self.config.scopes + ["global"]:
            for target_node, required_inputs in self.inputs[scope].items():
                compose = defaultdict(list)

                # Determine where each input comes from
                for req_input in required_inputs:
                    if self.outputs["global"].get(req_input):
                        source = self.outputs["global"][req_input][0]
                    elif self.outputs[scope].get(req_input):
                        source = self.outputs[scope][req_input][0]
                    elif req_input in self.nanoAOD_inputs:
                        source = f"global_{self.in_file_name}"
                    else:
                        raise ValueError(f"Input {req_input} is missing from NanoAOD and producers.")

                    compose[source].append(req_input)

                # Create edges grouped by source node
                for source_node, input_names in compose.items():
                    composite_name = ", ".join(sorted(input_names))
                    if self.debugging:
                        print(f"Adding Edge {composite_name} from {source_node} to {target_node}")
                    self.add_edge(source=source_node, target=target_node, scope=scope, name=composite_name)

    def parse_VectorProducer(self, producer, parent, scope, align=""):
        """Handles legacy vector producers by passing them to the regex parser."""
        if self.debugging:
            print(align + f"Adding VectorProducer: {producer.name}")
        print(f"!!! {producer.name} is a legacy producer and should be replaced with ExtendedVectorProducer !!!")

        if producer.call.startswith("event::filter::Flag"):
            self.parse_from_call(producer=producer, parent=parent, scope=scope, align=align)
        else:
            raise NotImplementedError(f"The call {producer.call} does not have legacy support.")

    def parse_ExtendedVectorProducer(self, producer, parent, scope, align=""):
        """Handles extended vector producers, mapping nested inputs and outputs."""
        if self.debugging:
            print(align + f"Adding ExtendedVectorProducer: {producer.name}")

        group_id = f"{producer.name}_v"
        self.add_node(id_name=group_id, name=producer.name, scope=scope, parent=parent, node_type="vector")
        vec_config = self.config.config_parameters[scope][producer.vec_config]

        for i_c, c in enumerate(vec_config):
            vector_id = f"{group_id}_{i_c}"
            vector_name = f"{producer.name}_{i_c}"

            if self.debugging:
                print(align + "    " + f"Adding Producer: {vector_name}")
            self.add_node(id_name=vector_id, name=vector_name, scope=scope, parent=group_id)

            if isinstance(producer.input[scope], list):
                for n in producer.input[scope]:
                    self.add_input(n.name, vector_id, scope)
                    if self.debugging:
                        print(align + "        " + f"Adding Input: {n.name}")

            if isinstance(producer.output, list):
                vec_output = c[producer.outputname]
                self.vec_output_mappings[producer.output_group.name] = producer.outputname
                self.add_output(vec_output, vector_id, scope)
                if self.debugging:
                    print(align + "        " + f"Adding Output: {vec_output}")

    def parse_ProducerGroup(self, producer, parent, scope, align=""):
        """Recursively parses a group of sub-producers, grouping them under a common parent."""
        if self.debugging:
            print(align + f"Adding ProducerGroup: {producer.name}")

        group_id = f"{producer.name}_g"
        self.add_node(id_name=group_id, name=producer.name, scope=scope, parent=parent)

        for p in producer.producers[scope]:
            self.parse_producer(producer=p, parent=group_id, scope=scope, align=align + "    ")

    def parse_Filter(self, producer, parent, scope, align=""):
        """Recursively parses a group of filters under a common parent group."""
        if self.debugging:
            print(align + f"Adding Filter Group: {producer.name}")

        group_id = f"{producer.name}_f"
        self.add_node(id_name=group_id, name=producer.name, scope=scope, parent=parent)

        for p in producer.producers[scope]:
            self.parse_producer(producer=p, parent=group_id, scope=scope, align=align + "    ")

    def parse_BaseFilter(self, producer, parent, scope, align=""):
        """Parses a base level filter node and records its required inputs."""
        if self.debugging:
            print(align + f"Adding BaseFilter: {producer.name}")
        print(f"!!! {producer.name} is a legacy producer and should be replaced with Filter !!!")

        self.add_node(id_name=producer.name, scope=scope, parent=parent)
        for n in producer.input[scope]:
            self.add_input(n.name, producer.name, scope)
            if self.debugging:
                print(align + "    " + f"Adding Input: {n.name}")

    def parse_Producer(self, producer, parent, scope, align=""):
        """Parses a standard producer node, mapping its simple list of inputs and outputs."""
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
        """
        Routing method. Determines the class type of the producer and dispatches
        it to the appropriate specific parsing method dynamically.
        """
        class_name = producer.__class__.__name__
        parser_method = getattr(self, f"parse_{class_name}", None)

        if parser_method:
            parser_method(producer, parent, scope, align)
        else:
            raise NotImplementedError(f"Unknown Producer class {class_name}")


def create_graph(configuration, nanoAOD_inputs, DAG_dir, json_name):
    """Entry point function to instantiate a GraphParser and execute generation."""
    graph = GraphParser(configuration, nanoAOD_inputs, DAG_dir, debugging=False)
    graph.generate_graph()
    for scope in graph.config.scopes:
        graph.save_graph(scope, json_name)
