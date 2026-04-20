from math import sqrt
import re
import json
import os
from collections import defaultdict
from code_generation.quantity import QuantityGroup
from code_generation.helpers import is_empty

def create_graph(configuration, nanoAOD_inputs, DAG_dir, json_name, debugging=False):
    """Entry point function to instantiate a GraphParser and execute generation."""

    for active_scope in configuration.scopes:
        if active_scope == "global":
            pass
        else:
            if hasattr(configuration, "input_quantities_mapping"):
                CROWN_input_quantities = configuration.input_quantities_mapping[
                    active_scope
                ][""]
            else:
                CROWN_input_quantities = None
            graph = GraphParser(
                configuration,
                nanoAOD_inputs,
                active_scope,
                DAG_dir,
                CROWN_input_quantities=CROWN_input_quantities,
                debugging=debugging,
            )
            graph.generate_graph()
            graph.save_graph(json_name)


class GraphParser:
    """
    Parses configuration data to generate Directed Acyclic Graphs (DAGs)
    representing the dependencies and flow of producers, filters, and I/O.
    """

    def __init__(
        self,
        config,
        nanoAOD_inputs,
        active_scope,
        DAG_dir=None,
        CROWN_input_quantities=None,
        debugging=False,
    ):
        self.config = config
        self.active_scope = active_scope
        if active_scope == "global":
            self.scopes = ["global"]
        else:
            self.scopes = ["global", self.active_scope]
        self.nanoAOD_inputs = nanoAOD_inputs
        if CROWN_input_quantities is None:
            self.CROWN_input_quantities = {}
        else:
            self.CROWN_input_quantities = CROWN_input_quantities
        self.debugging = debugging
        self.nodes = defaultdict()
        self.connections = defaultdict(lambda: defaultdict(list))
        self.edges = []
        self.inputs = defaultdict(lambda: defaultdict(list))
        self.outputs = defaultdict(lambda: defaultdict(list))
        self.vec_output_mappings = {}
        self.node_family_register = defaultdict(
            lambda: {
                "file_in": {"NanoAOD": [], "Ntuple": []},
                "file_out": [],
                "node_call_data": {},
            }
        )
        self.edge_family_register = defaultdict(
            lambda: {"members": [], "familyHead": None}
        )
        self.DAG_data = None
        self.DAG_dir = DAG_dir

    def generate_graph(self):
        """
        The main orchestrator. Iterates through all scopes of the configuration,
        generates the nodes and edges, verifies output ambiguity, and generates metadata.
        """
        for scope in self.scopes:
            if not is_empty(self.config.producers[scope]):
                if self.debugging:
                    print(f"\nFor scope {scope}:")
                self.add_node(
                    id_name="scope",
                    name=f"{scope} scope",
                    scope=scope,
                    node_type="scope",
                )

            # Parse all producers in this scope
            for p in self.config.producers[scope]:
                self.parse_Producer_routing(
                    producer=p, parent="scope", scope=scope, align="    "
                )

            # Check outputs for ambiguous assignments (multiple origins for the same output)
            all_keys = list(
                set(
                    list(self.outputs["global"].keys())
                    + list(self.outputs[scope].keys())
                )
            )
            for key in all_keys:
                total_output_nodes = self.outputs["global"].get(key, []) + self.outputs[
                    scope
                ].get(key, [])
                if len(set(total_output_nodes)) != 1:
                    raise ValueError(
                        f"Output {key} has multiple origins: {total_output_nodes}"
                    )

            # Determine Output nodes by tracing configured scope outputs back to their producers
            if hasattr(self.config, "outputs") and scope in self.config.outputs:
                for output in self.config.outputs[scope]:
                    if isinstance(output, QuantityGroup):
                        if output.vec_config in self.config.config_parameters.get(
                            scope, {}
                        ):
                            vec_config = self.config.config_parameters[scope][
                                output.vec_config
                            ]
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

        # Identify nodes without a specific type and without any outgoing connections
        # This excludes filters, groups, vector groups, and scopes.
        # While they are only defined afterwards this also excludes proxy nodes and edges (branch, twig, leaf).
        for node_id, node_data in self.nodes.items():
            if (
                not node_data.get("type")
                and not node_data.get("is_out")
                and node_id not in self.connections
            ):
                node_data["type"] = "stump"

        # Bundle edges for quantities with multiple targets
        self.bundle_edges()

        # Prepare node formatting for usage with cytoscape
        formatted_nodes = [{"data": n} for n in self.nodes.values()]

        # Compile final DAG data
        # self.DAG_data = formatted_nodes + self.edges
        self.DAG_data = sorted(
            formatted_nodes, key=lambda d: d["data"].get("label", "")
        ) + sorted(self.edges, key=lambda d: d["data"]["target"])

    def parse_Producer_routing(self, producer, parent, scope, align=""):
        """
        Router to producer classes.
        """
        class_name = producer.__class__.__name__

        if class_name == "VectorProducer":
            self.parse_VectorProducer(
                producer=producer, parent=parent, scope=scope, align=align
            )
        elif class_name == "ExtendedVectorProducer":
            self.parse_ExtendedVectorProducer(
                producer=producer, parent=parent, scope=scope, align=align
            )
        elif class_name == "ProducerGroup":
            self.parse_ProducerGroup(
                producer=producer, parent=parent, scope=scope, align=align
            )
        elif class_name == "Filter":
            self.parse_Filter(
                producer=producer, parent=parent, scope=scope, align=align
            )
        elif class_name == "BaseFilter":
            self.parse_BaseFilter(
                producer=producer, parent=parent, scope=scope, align=align
            )
        elif class_name == "Producer":
            self.parse_Producer(
                producer=producer, parent=parent, scope=scope, align=align
            )
        else:
            raise NotImplementedError(f"Unknown Producer class {class_name}")

    def parse_VectorProducer(self, producer, parent, scope, align=""):
        """
        Uses regex to extract vector configurations from a legacy producer's string call
        and maps them into graph nodes and inputs.
        Currently only supports 'event::filter::Flag'.
        """
        if self.debugging:
            print(align + f"Adding VectorProducer: {producer.name}")
        print(
            f"!!! {producer.name} is a legacy producer and should be replaced with ExtendedVectorProducer !!!"
        )

        if producer.call.startswith("event::filter::Flag"):
            self.parse_Flag_from_call(
                producer=producer, parent=parent, scope=scope, align=align
            )
        else:
            raise NotImplementedError(
                f"The call {producer.call} does not have legacy support."
            )

    def parse_Flag_from_call(self, producer, parent, scope, align=""):
        """
        Uses regex to extract vector configurations from a legacy producer's string call
        and maps them onto configs.
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
            raise ValueError(
                f"Input name from {call} not in producer vector configs {producer.vec_configs}"
            )

        vec_input_index = producer.vec_configs.index(input_vec_config)
        call = producer.call
        self.add_node(
            id_name=group_id,
            name=producer.name,
            scope=scope,
            parent=parent,
            node_type="vector",
        )
        call = producer.call
        vec_configs = [
            self.config.config_parameters[scope][c] for c in producer.vec_configs
        ]

        for i_c, c in enumerate(zip(*vec_configs)):
            input_name = c[vec_input_index]
            vector_id = f"{group_id}_{i_c}"
            vector_name = f"{producer.name}_{i_c}"
            vector_config_dict = {vc: cc for vc, cc in zip(producer.vec_configs, c)}
            config_data = self.extract_configs(call, scope, vector_config_dict)
            if self.debugging:
                print(align + "    " + f"Adding Filter: {vector_name}")
            self.add_node(
                id_name=vector_id,
                name=vector_name,
                scope=scope,
                parent=group_id,
                is_filter=producer.is_filter,
                node_call_data={"call": call, "configs": config_data},
            )
            self.add_input(input_name, vector_id, scope)

            if isinstance(producer.output, list):
                raise NotImplementedError(
                    "List outputs for legacy parsed calls are not supported."
                )

    def parse_ExtendedVectorProducer(self, producer, parent, scope, align=""):
        """
        Parses the ExtendedVectorProducer class into graph nodes and inputs/outputs
        """
        if self.debugging:
            print(align + f"Adding ExtendedVectorProducer: {producer.name}")

        group_id = f"{producer.name}_v"
        self.add_node(
            id_name=group_id,
            name=producer.name,
            scope=scope,
            parent=parent,
            node_type="vector",
        )
        call = producer.call
        vec_config = self.config.config_parameters[scope][producer.vec_config]

        for i_c, c in enumerate(vec_config):
            vector_id = f"{group_id}_{i_c}"
            vector_name = f"{producer.name}_{i_c}"
            config_data = self.extract_configs(call, scope, c)
            if self.debugging:
                print(align + "    " + f"Adding Producer: {vector_name}")
            self.add_node(
                id_name=vector_id,
                name=vector_name,
                scope=scope,
                parent=group_id,
                is_filter=producer.is_filter,
                node_call_data={"call": call, "configs": config_data},
            )

            if isinstance(producer.input[scope], list):
                for n in set(producer.input[scope]):
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
        """
        Parses the ProducerGroup class into graph nodes and inputs/outputs.
        Calls parse_Producer_routing for each producer in the group.
        """
        if self.debugging:
            print(align + f"Adding ProducerGroup: {producer.name}")

        group_id = f"{producer.name}_g"
        self.add_node(
            id_name=group_id,
            name=producer.name,
            scope=scope,
            parent=parent,
            node_type="group",
        )

        for p in producer.producers[scope]:
            self.parse_Producer_routing(
                producer=p, parent=group_id, scope=scope, align=align + "    "
            )

    def parse_Filter(self, producer, parent, scope, align=""):
        """
        Parses the Filter class into graph nodes and inputs.
        Calls parse_Producer_routing for each producer in the group,
        similar to parse_ProducerGroup.
        """
        if self.debugging:
            print(align + f"Adding Filter Group: {producer.name}")

        group_id = f"{producer.name}_f"
        self.add_node(
            id_name=group_id,
            name=producer.name,
            scope=scope,
            parent=parent,
            node_type="group",
        )

        for p in producer.producers[scope]:
            self.parse_Producer_routing(
                producer=p, parent=group_id, scope=scope, align=align + "    "
            )

    def parse_BaseFilter(self, producer, parent, scope, align=""):
        """
        Parses the Legacy BaseFilter class into graph nodes and inputs.
        """
        if self.debugging:
            print(align + f"Adding BaseFilter: {producer.name}")
        print(
            f"!!! {producer.name} is a legacy producer and should be replaced with Filter !!!"
        )
        call = producer.call
        config_data = self.extract_configs(call, scope)
        self.add_node(
            id_name=producer.name,
            scope=scope,
            parent=parent,
            is_filter=producer.is_filter,
            node_call_data={"call": call, "configs": config_data},
        )
        for n in set(producer.input[scope]):
            self.add_input(n.name, producer.name, scope)
            if self.debugging:
                print(align + "    " + f"Adding Input: {n.name}")

    def parse_Producer(self, producer, parent, scope, align=""):
        """
        Parses the Producer class into graph nodes and inputs/outputs.
        """
        if self.debugging:
            print(align + f"Adding Producer: {producer.name}")
        call = producer.call
        config_data = self.extract_configs(call, scope)
        self.add_node(
            id_name=producer.name,
            scope=scope,
            parent=parent,
            is_filter=producer.is_filter,
            node_call_data={"call": call, "configs": config_data},
        )

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

    def set_is_out(self, scope, req_out):
        """Sets the node as a source of an Ntuple output."""
        producers = self.outputs[scope].get(req_out, []) + self.outputs["global"].get(
            req_out, []
        )
        if len(producers) > 0:
            if len(producers) != 1:
                raise ValueError(f"Num producers for out {req_out}: {len(producers)}")
            if self.nodes.get(producers[0]):
                self.nodes[producers[0]]["is_out"] = True
                self.node_family_register[producers[0]]["file_out"].append(req_out)
            else:
                raise ValueError(
                    f"Source {producers[0]} is neither part of {scope} nor global scope."
                )
        else:
            if not (
                req_out in self.nanoAOD_inputs or req_out in self.CROWN_input_quantities
            ):
                raise ValueError(
                    f"Requested output {req_out} is neither provided by a producer nor NanoAOD/Ntuple."
                )

    def assemble_connections(self):
        """
        Resolves the inputs vs. outputs mappings to draw the actual connecting
        lines (edges) between nodes in the DAG.
        Determines the source of each input: Producer, NanoAOD, or Ntuple.
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
                        self.nodes[target_node].setdefault("is_in", set()).add(
                            "NanoAOD"
                        )
                        self.node_family_register[target_node]["file_in"][
                            "NanoAOD"
                        ].append(req_input)
                    elif req_input in self.CROWN_input_quantities:
                        self.nodes[target_node].setdefault("is_in", set()).add("Ntuple")
                        self.node_family_register[target_node]["file_in"][
                            "Ntuple"
                        ].append(req_input)
                    else:
                        raise ValueError(
                            f"Input {req_input} is missing from NanoAOD and producers."
                        )

                # Create edges grouped by source node
                for source_node, input_names in compose.items():
                    for input_name in input_names:
                        if self.debugging:
                            print(
                                f"Adding Connection {input_name} from {source_node} to {target_node}"
                            )
                        self.add_connection(
                            source=source_node, target=target_node, name=input_name
                        )
        for node in self.nodes.values():
            if "is_in" in node:
                node["is_in"] = ",".join(sorted(list(node["is_in"])))

    def bundle_edges(self):
        """
        Bundles connections based on common target location.
        Creates proxy nodes and edges with width based on number of connections.
        """
        for source, dt in self.connections.items():
            for name, targets in dt.items():
                if len(targets) == 1:
                    target = targets[0]
                    # edge_id = f"{name}_{target}"
                    location = self.nodes[source].get("parent")
                    proxy_node_name = f"proxy_{name}_{location}"
                    proxy_edge_name = f"proxyedge_{name}_{location}"
                    leaf_edge_name = f"{name}_{target}"
                    if self.debugging:
                        print(f"Label at {location}")
                        print(
                            f"Adding proxy node {proxy_node_name} with parent {location}"
                        )
                        print(
                            f"Adding proxy edge {proxy_edge_name} from {source} to {proxy_node_name} with weight 1"
                        )
                        print(
                            f"Adding final edge {leaf_edge_name} from {proxy_node_name} to {target} with weight 1"
                        )
                    self.add_node(
                        id_name=proxy_node_name,
                        name=name,
                        parent=location,
                        node_type="proxy",
                        family=f"edge_{name}",
                    )
                    self.add_edge(
                        source=source,
                        target=proxy_node_name,
                        edge_id=proxy_edge_name,
                        name=name,
                        weight=1,
                        edge_type="twig",
                        family=f"edge_{name}",
                    )
                    self.add_edge(
                        source=proxy_node_name,
                        target=target,
                        edge_id=leaf_edge_name,
                        name=name,
                        weight=1,
                        edge_type="leaf",
                        family=f"edge_{name}",
                    )
                else:
                    source_history = self.get_ancestors(source)
                    relative_targets = [
                        self.get_relative(source_history, self.get_ancestors(target))
                        for target in targets
                    ]
                    start_idx = list(range(len(relative_targets)))
                    self.construct_proxies(
                        name,
                        relative_targets,
                        start_idx,
                        start_idx,
                        0,
                        relative_targets[0][0],
                    )

    def get_ancestors(self, node_name):
        """Returns recursive list of parents of a given node."""
        parent = self.nodes[node_name].get("parent")
        if not parent:
            return ["root", node_name]
        else:
            return self.get_ancestors(parent) + [node_name]

    def get_relative(self, source_hist, target_hist):
        """Returns the relative path from source to target."""
        min_len = min(len(source_hist), len(target_hist))
        for i in range(min_len):
            if source_hist[i] != target_hist[i]:
                return (
                    list(reversed(source_hist[i:]))
                    + [source_hist[i - 1]]
                    + target_hist[i:]
                )

    def construct_proxies(self, name, data, valid_idx, tot_idx, rank, source):
        """Recursive function to create proxy nodes and edges."""
        data_slice = [data[i][rank] for i in tot_idx if i in valid_idx]
        if len(set(data_slice)) <= 1:
            self.construct_proxies(name, data, valid_idx, tot_idx, rank + 1, source)
        else:
            proxy_loc = data[valid_idx[0]][rank - 1]
            proxy_node_name = f"proxy_{name}_{proxy_loc}"
            proxy_edge_name = f"proxyedge_{name}_{proxy_loc}"
            if self.debugging:
                print(f"Junktion at {proxy_loc}")
                print(f"Adding proxy node {proxy_node_name} with parent {proxy_loc}")
                print(
                    f"Adding proxy edge {proxy_edge_name} from {source} to {proxy_node_name} with weight {len(valid_idx)}"
                )
            self.add_node(
                id_name=proxy_node_name,
                name=name,
                parent=proxy_loc,
                node_type="proxy",
                family=f"edge_{name}",
            )
            self.add_edge(
                source=source,
                target=proxy_node_name,
                edge_id=proxy_edge_name,
                name=name,
                weight=len(valid_idx),
                edge_type="branch",
                family=f"edge_{name}",
            )
            for val in set(data_slice):
                part_idx = [valid_idx[i] for i, d in enumerate(data_slice) if d == val]
                if len(part_idx) <= 1:
                    leaf_name = data[part_idx[0]][-1]
                    leaf_edge_name = f"{name}_{leaf_name}"
                    if self.debugging:
                        print(
                            f"Adding final edge {leaf_edge_name} from {proxy_node_name} to {leaf_name} with weight 1"
                        )
                    self.add_edge(
                        source=proxy_node_name,
                        target=leaf_name,
                        edge_id=leaf_edge_name,
                        name=name,
                        weight=1,
                        edge_type="leaf",
                        family=f"edge_{name}",
                    )
                else:
                    self.construct_proxies(
                        name, data, part_idx, tot_idx, rank + 1, proxy_node_name
                    )

    def extract_configs(self, call, scope, vector_configs=None, ignore={}):
        """Extracts configuration parameters from a call string."""
        if is_empty(ignore):
            ignore = {
                "df",
                "output",
                "input",
                "output_vec",
                "input_vec",
                "vec_open",
                "vec_close",
            }
        else:
            ignore = set(ignore)

        pattern = r"\{(\w+)\}"
        matches = re.findall(pattern, call)

        config_parameters = [m for m in matches if m not in ignore]

        config_dict = {}
        for c in config_parameters:
            if vector_configs != None and not is_empty(vector_configs.get(c)):
                config_dict[c] = vector_configs[c]
            elif not is_empty(self.config.config_parameters[scope].get(c)):
                config_dict[c] = self.config.config_parameters[scope][c]
            else:
                raise ValueError(f"Unknown config parameter {c}")
        return config_dict

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

    def add_node(
        self,
        id_name,
        name=None,
        scope=None,
        parent=None,
        node_type=None,
        is_filter=False,
        family=None,
        node_call_data=None,
    ):
        """
        Creates a standardized node object and appends it to the graph's node list.
        Also updates the node family register and edge family register.
        """

        if node_type and is_filter:
            raise ValueError("Cannot specify both node_type and is_filter")
        if is_filter:
            node_type = "filter"

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
        # Any node with a family is assumed to be a proxy node part of an edge connection
        if family:
            add_data["family"] = family
            self.edge_family_register[family]["members"].append(full_id)
            if not self.edge_family_register[family]["familyHead"]:
                self.edge_family_register[family]["familyHead"] = full_id
        else:
            # Standalone nodes will only ever have one member in their family, themselves
            if self.node_family_register.get(full_id):
                raise ValueError(f"Node {full_id} already exists in family register")
            add_data["family"] = full_id
            if node_call_data:
                self.node_family_register[full_id]["node_call_data"] = node_call_data

        self.nodes[full_id] = add_data

    def add_connection(self, source, target, name):
        """Creates a directional edge between a source node and a target node."""
        self.connections[source][name].append(target)

    def add_edge(self, source, target, edge_id, name, weight, edge_type, family):
        """Creates a directional edge between a source node and a target node."""
        self.edge_family_register[family]["members"].append(edge_id)
        self.edges.append(
            {
                "data": {
                    "id": edge_id,
                    "source": source,
                    "target": target,
                    "type": edge_type,
                    "width": sqrt(weight),
                    "description": name,
                    "family": family,
                }
            }
        )

    def update_DAG_file_list(self, config_path, new_era, new_sample, new_scope):
        """
        Updates the master JSON tracking file to ensure new eras, samples,
        and scopes are registered without duplicates.
        """
        if os.path.exists(config_path):
            with open(config_path, "r") as f:
                try:
                    data = json.load(f)
                except json.JSONDecodeError:
                    data = {"era": [], "sample": [], "scope": []}
        else:
            data = {"era": [], "sample": [], "scope": []}

        data["era"] = sorted(list(set(data["era"] + [str(new_era)])))
        data["sample"] = sorted(list(set(data["sample"] + [str(new_sample)])))
        data["scope"] = sorted(list(set(data["scope"] + [str(new_scope)])))

        with open(config_path, "w") as f:
            json.dump(data, f, indent=4)

        if self.debugging:
            print(f"Updated {config_path} with {new_era}/{new_sample}/{new_scope}")

    def save_graph(self, name):
        """
        Compiles the global and scoped DAG data and saves it to a JSON file.
        Also triggers an update to the master DAG file list.
        """
        path = f"{name}_{self.config.era}_{self.config.sample}_{self.active_scope}.json"
        if self.DAG_dir:
            path = os.path.join(self.DAG_dir, path)
            os.makedirs(self.DAG_dir, exist_ok=True)

        meta_data = {
            "nodeFamilyRegister": self.node_family_register,
            "edgeFamilyRegister": self.edge_family_register,
        }
        full_dict = {"elementData": self.DAG_data, "metaData": meta_data}

        with open(path, "w") as f:
            json.dump(full_dict, f, indent=4)  # , sort_keys=True)

        self.update_DAG_file_list(
            os.path.join(self.DAG_dir, "DAG_files.json"),
            self.config.era,
            self.config.sample,
            self.active_scope,
        )

