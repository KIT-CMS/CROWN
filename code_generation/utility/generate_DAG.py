from math import sqrt
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
    def __init__(self, config, nanoAOD_inputs, active_scope, DAG_dir=None, debugging=False):
        self.config = config
        self.active_scope = active_scope
        if active_scope == "global":
            self.scopes = ["global"]
        else:
            self.scopes = ["global", self.active_scope]
        self.nanoAOD_inputs = nanoAOD_inputs
        self.debugging = debugging
        self.nodes = defaultdict()
        self.connections = defaultdict(
            lambda: defaultdict(list)
        )
        self.edges = []
        self.inputs = defaultdict(lambda: defaultdict(list))
        self.outputs = defaultdict(lambda: defaultdict(list))
        self.vec_output_mappings = {}

        self.DAG_data = {}
        self.DAG_dir = DAG_dir

    def save_graph(self, name):
        """
        Compiles the global and scoped DAG data and saves it to a JSON file.
        Also triggers an update to the master DAG file list.
        """
        path = f"{name}_{self.config.era}_{self.config.sample}_{self.active_scope}.json"
        if self.DAG_dir:
            path = os.path.join(self.DAG_dir, path)
            os.makedirs(self.DAG_dir, exist_ok=True)

        with open(path, "w") as f:
            json.dump(self.DAG_data, f, indent=4, sort_keys=True)

        self.update_DAG_file_list(
            os.path.join(self.DAG_dir, "DAG_files.json"),
            self.config.era,
            self.config.sample,
            self.active_scope
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

        if self.debugging:
            print(f"Updated {config_path} with {new_era}/{new_sample}/{new_scope}")

    def generate_graph(self):
        """
        The main orchestrator. Iterates through all scopes of the configuration,
        generates the inputs/outputs/nodes, verifies output ambiguity, and creates edges.
        """
        for scope in self.scopes:
            if self.debugging:
                print(f"\nFor scope {scope}:")
            self.add_node(id_name="scope", name=f"{scope} scope", scope=scope, node_type="scope")

            # Parse all producers in this scope
            for p in self.config.producers[scope]:
                self.parse_producer(producer=p, parent="scope", scope=scope, align="    ")

            # Check outputs for ambiguous assignments (multiple origins for the same output)
            all_keys = list(set(list(self.outputs["global"].keys()) + list(self.outputs[scope].keys())))
            for key in all_keys:
                total_output_nodes = self.outputs["global"].get(key, []) + self.outputs[scope].get(key, [])
                if len(set(total_output_nodes)) != 1:
                    raise ValueError(f"Output {key} has multiple origins: {total_output_nodes}")

            # Determine Output nodes by tracing configured scope outputs back to their producers
            if hasattr(self.config, 'outputs') and scope in self.config.outputs:
                for output in self.config.outputs[scope]:
                    if isinstance(output, QuantityGroup):
                        if output.vec_config in self.config.config_parameters.get(scope, {}):
                            vec_config = self.config.config_parameters[scope][output.vec_config]
                            vec_output_name = self.vec_output_mappings.get(output.name)
                            if vec_output_name:
                                for o in vec_config:
                                    req_out = o.get(vec_output_name)
                                    if req_out:
                                        self.set_is_out(scope, req_out)
                    else:
                        self.set_is_out(scope, output.name)
        # Derive connections from inputs/outputs
        self.assemble_connections()

        # Bundle edges for quantities with multiple targets
        self.bundle_edges()

        #Prepare node formatting for usage with cytoscape
        formatted_nodes = [{"data": n} for n in self.nodes.values()]

        # Compile final DAG data
        self.DAG_data = formatted_nodes + self.edges

    def set_is_out(self, scope, req_out):
        producers = self.outputs[scope].get(req_out, []) + self.outputs["global"].get(req_out, [])
        if len(producers)>0:
            assert len(producers)==1, f"Num producers for out {req_out}: {len(producers)}"
            if self.nodes.get(producers[0]):
                self.nodes[producers[0]]["is_out"] = True
            else:
                raise ValueError(f"Source {producers[0]} is neither part of {scope} nor global scope.")
        else:
            if not req_out in self.nanoAOD_inputs:
                raise ValueError(f"Requested output {req_out} is neither provided by a producer nor NanoAOD.")

    def assemble_connections(self):
        """
        Resolves the inputs vs. outputs mappings to draw the actual connecting
        lines (edges) between nodes in the DAG. Records IN scope nodes.
        """
        for scope in self.scopes:
            for target_node, required_inputs in self.inputs[scope].items():
                compose = defaultdict(list)

                # Determine where each input comes from
                for req_input in required_inputs:
                    if self.outputs["global"].get(req_input):
                        source = self.outputs["global"][req_input][0]
                        compose[source].append(req_input)
                    elif self.outputs[scope].get(req_input):
                        source = self.outputs[scope][req_input][0]
                        compose[source].append(req_input)
                    elif req_input in self.nanoAOD_inputs:
                        self.nodes[target_node]["is_in"] = True
                    else:
                        raise ValueError(f"Input {req_input} is missing from NanoAOD and producers.")

                # Create edges grouped by source node
                for source_node, input_names in compose.items():
                    for input_name in input_names:
                        if self.debugging:
                            print(f"Adding Connection {input_name} from {source_node} to {target_node}")
                        self.add_connection(source=source_node, target=target_node, name=input_name)

    def bundle_edges(self):
        # self.connections[scope][name][source].append(target)
        for source, dt in self.connections.items():
            for name, targets in dt.items():
                if len(targets)==1:
                    target = targets[0]
                    edge_id = f"{name}_{target}"
                    self.add_edge(source, targets[0], edge_id, name, 1, "leaf")
                else:
                    source_history = self.get_ancestors(source)
                    relative_targets = [self.get_relative(source_history, self.get_ancestors(target)) for target in targets]
                    start_idx = list(range(len(relative_targets)))
                    self.construct_proxies(name, relative_targets, start_idx, start_idx, 0, relative_targets[0][0])

    def get_ancestors(self, node_name):
        parent = self.nodes[node_name].get("parent")
        if not parent:
            return ["root", node_name]
        else:
            return self.get_ancestors(parent) + [node_name]

    def get_relative(self, source_hist, target_hist):
        min_len = min(len(source_hist), len(target_hist))
        for i in range(min_len):
            if source_hist[i] != target_hist[i]:
                return list(reversed(source_hist[i:])) + [source_hist[i-1]] + target_hist[i:]


    def construct_proxies(self, name, data, valid_idx, tot_idx, rank, source):
        data_slice = [data[i][rank] for i in tot_idx if i in valid_idx]
        if len(set(data_slice)) <= 1:
            self.construct_proxies(name, data, valid_idx, tot_idx, rank + 1, source)
        else:
            proxy_loc = data[valid_idx[0]][rank-1]
            proxy_node_name = f"proxy_{name}_{proxy_loc}"
            proxy_edge_name = f"proxyedge_{name}_{proxy_loc}"
            if self.debugging:
                print(f"Junktion at {proxy_loc}")
                print(f"Adding proxy node {proxy_node_name} with parent {proxy_loc}")
                print(f"Adding proxy edge {proxy_edge_name} from {source} to {proxy_node_name} with weight {len(valid_idx)}")
            self.add_node(proxy_node_name, parent=proxy_loc, node_type="proxy")
            self.add_edge(source, proxy_node_name, proxy_edge_name, name, len(valid_idx), "branch")
            for val in set(data_slice):
                part_idx = [valid_idx[i] for i, d in enumerate(data_slice) if d == val]
                if len(part_idx) <= 1:
                    leaf_name = data[part_idx[0]][-1]
                    leaf_edge_name = f"{name}_{leaf_name}"
                    if self.debugging:
                        print(f"Adding final edge {leaf_edge_name} from {proxy_node_name} to {leaf_name} with weight 1")
                    self.add_edge(proxy_node_name, leaf_name, leaf_edge_name, name, 1, "leaf")
                else:
                    self.construct_proxies(name, data, part_idx, tot_idx, rank + 1, proxy_node_name)

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

    def parse_VectorProducer(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding VectorProducer: {producer.name}")
        print(f"!!! {producer.name} is a legacy producer and should be replaced with ExtendedVectorProducer !!!")

        if producer.call.startswith("event::filter::Flag"):
            self.parse_from_call(producer=producer, parent=parent, scope=scope, align=align)
        else:
            raise NotImplementedError(f"The call {producer.call} does not have legacy support.")

    def parse_ExtendedVectorProducer(self, producer, parent, scope, align=""):
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
                for n in set(producer.input[scope]):
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
        if self.debugging:
            print(align + f"Adding ProducerGroup: {producer.name}")

        group_id = f"{producer.name}_g"
        self.add_node(id_name=group_id, name=producer.name, scope=scope, parent=parent)

        for p in producer.producers[scope]:
            self.parse_producer(producer=p, parent=group_id, scope=scope, align=align + "    ")

    def parse_Filter(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding Filter Group: {producer.name}")

        group_id = f"{producer.name}_f"
        self.add_node(id_name=group_id, name=producer.name, scope=scope, parent=parent)

        for p in producer.producers[scope]:
            self.parse_producer(producer=p, parent=group_id, scope=scope, align=align + "    ")

    def parse_BaseFilter(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding BaseFilter: {producer.name}")
        print(f"!!! {producer.name} is a legacy producer and should be replaced with Filter !!!")

        self.add_node(id_name=producer.name, scope=scope, parent=parent)
        for n in set(producer.input[scope]):
            self.add_input(n.name, producer.name, scope)
            if self.debugging:
                print(align + "    " + f"Adding Input: {n.name}")

    def parse_Producer(self, producer, parent, scope, align=""):
        if self.debugging:
            print(align + f"Adding Producer: {producer.name}")

        self.add_node(id_name=producer.name, scope=scope, parent=parent)

        if isinstance(producer.input[scope], list):
            for n in set(producer.input[scope]):
                self.add_input(n.name, producer.name, scope)
                if self.debugging:
                    print(align + "    " + f"Adding Input: {n.name}")

        if isinstance(producer.output, list):
            for n in set(producer.output):
                self.add_output(n.name, producer.name, scope)
                if self.debugging:
                    print(align + "    " + f"Adding Output: {n.name}")

    def parse_producer(self, producer, parent, scope, align=""):
        class_name = producer.__class__.__name__
        parser_method = getattr(self, f"parse_{class_name}", None)

        if parser_method:
            parser_method(producer, parent, scope, align)
        else:
            raise NotImplementedError(f"Unknown Producer class {class_name}")

    def add_output(self, output_name, output_node, scope):
        """Registers a node as the creator of a specific output."""
        node_id = f"{scope}_{output_node}"
        if node_id in self.outputs[scope][output_name]:
            raise ValueError(f"Node {node_id} already exists in output nodes.")
        self.outputs[scope][output_name].append(node_id)

    def add_input(self, input_name, input_node, scope):
        """Registers a node's requirement for a specific input."""
        scoped_input = f"{scope}_{input_node}"
        if input_name in self.inputs[scope][scoped_input]:
            raise ValueError(f"Input {input_name} already exists in inputs.")
        self.inputs[scope][scoped_input].append(input_name)

    def add_node(self, id_name, name=None, scope=None, parent=None, node_type=None):
        """Creates a standardized node object and appends it to the graph's node list."""
        if not name:
            name = id_name
        if scope:
            full_id = f"{scope}_{id_name}"
        else:
            full_id = id_name
        add_data = {"id": full_id, "label": name}
        if parent:
            if scope:
                add_data["parent"] = f"{scope}_{parent}"
            else:
                add_data["parent"] = parent
        if node_type:
            add_data["type"] = node_type

        self.nodes[full_id] = add_data

    def add_connection(self, source, target, name):
        """Creates a directional edge between a source node and a target node."""
        self.connections[source][name].append(target)

    def add_edge(self, source, target, edge_id, name, weight, edge_type, family=None):
        """Creates a directional edge between a source node and a target node."""
        if not family:
            family = [edge_id]
        self.edges.append({
            "data": {
                "id":edge_id,
                "source": source,
                "target": target,
                "type": edge_type,
                "width": sqrt(weight),
                "description": name,
                "family": family
            }
        })

def create_graph(configuration, nanoAOD_inputs, DAG_dir, json_name, debugging=False):
    """Entry point function to instantiate a GraphParser and execute generation."""
    for active_scope in configuration.scopes:
        if active_scope=="global":
            pass
        else:
            graph = GraphParser(configuration, nanoAOD_inputs, active_scope, DAG_dir, debugging=debugging)
            graph.generate_graph()
            graph.save_graph(json_name)
