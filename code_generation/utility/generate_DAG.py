from math import sqrt
import re
import json
import os
from collections import defaultdict
from code_generation.quantity import QuantityGroup
from code_generation.helpers import is_empty


def create_graph(configuration, external_inputs, DAG_dir, json_name, debugging=False):
    """Instantiate a GraphParser and execute DAG generation.

    Args:
        configuration (object): The configuration object containing scopes and producers.
        external_inputs (list): List of external input quantities from NanoAOD.
        DAG_dir (str): Directory where the generated DAG JSON files will be saved.
        json_name (str): The base name for the output JSON file.
        debugging (bool, optional): If True, prints detailed debugging information. Defaults to False.
    """
    for active_scope in configuration.scopes:
        # Skip global scope as it is included in all other scopes
        if active_scope == "global":
            pass
        else:
            # Add Ntuple quantities for friends
            is_friend_config = False
            if hasattr(configuration, "input_quantities_mapping"):
                external_inputs = configuration.input_quantities_mapping[
                    active_scope
                ][""]
                is_friend_config = True
            graph = GraphParser(
                configuration,
                external_inputs,
                active_scope,
                DAG_dir,
                is_friend_config=is_friend_config,
                debugging=debugging,
            )
            # Generate and save graph for each scope separately
            graph.generate_graph()
            graph.save_graph(json_name)


class GraphParser:
    """Parses configuration data to generate Directed Acyclic Graphs (DAGs).

    Represents the dependencies and flow of producers, filters, and I/O.
    """

    def __init__(
        self,
        config,
        external_inputs,
        active_scope,
        DAG_dir="",
        is_friend_config=False,
        debugging=False,
    ):
        """Initialize the GraphParser object.

        Args:
            config (object): The configuration object.
            external_inputs (list): List of external inputs required by the graph.
            active_scope (str): The specific scope being parsed in addition to the global scope.
            DAG_dir (str, optional): The directory for saving DAG files. Defaults to "".
            is_friend_config (bool, optional): Indicates if it's a friend configuration. Defaults to False.
            debugging (bool, optional): Enables debug print statements. Defaults to False.
        """
        self.config = config
        self.active_scope = active_scope
        if active_scope == "global":
            self.scopes = ["global"]
        else:
            self.scopes = ["global", self.active_scope]
        self.external_inputs = external_inputs
        self.is_friend_config = is_friend_config
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
        """Execute the main orchestration to generate nodes and edges.

        Iterates through all scopes of the configuration, builds the DAG components,
        verifies output ambiguity, bundles edges, and generates formatting metadata
        for visualization tools like Cytoscape.
        """
        # Determine producers from configuration for global and active scope
        for scope in self.scopes:
            self.parse_scope_producers(scope)

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

        # Compile final DAG data and sort by label
        self.DAG_data = sorted(
            formatted_nodes, key=lambda d: d["data"].get("label", "")
        ) + sorted(self.edges, key=lambda d: d["data"]["target"])

    def parse_scope_producers(self, scope):
        """Parse producers and outputs for a specific scope.

        Args:
            scope (str): The scope name to process.

        Raises:
            ValueError: If an output has multiple origins (ambiguous assignments).
        """
        if not is_empty(self.config.producers[scope]):
            if self.debugging:
                print(f"\nFor scope {scope}:")
            # Add top level scope node
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

        # Determine nodes writing to Ntuple by tracing configured scope outputs back to their producers
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

    def parse_Producer_routing(self, producer, parent, scope, align=""):
        """Route parsing logic based on the specific producer class type.

        Args:
            producer (object): The producer instance to be parsed.
            parent (str): The ID of the parent node.
            scope (str): The active scope.
            align (str, optional): Spacing string used for debug print alignment. Defaults to "".

        Raises:
            NotImplementedError: If the producer's class is unknown.
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
        """Map legacy vector configurations to graph nodes and inputs using regex.

        Note:
            Currently only supports 'event::filter::Flag'. This is a legacy function
            and producers should be migrated to ExtendedVectorProducer.

        Args:
            producer (object): The legacy vector producer instance.
            parent (str): The ID of the parent node.
            scope (str): The active scope.
            align (str, optional): Debug print alignment spacing. Defaults to "".

        Raises:
            NotImplementedError: If the producer call does not have legacy support.
        """
        if self.debugging:
            print(align + f"Adding VectorProducer: {producer.name}")
        print(
            f"!!! {producer.name} is a legacy producer and should be replaced with ExtendedVectorProducer !!!"
        )

        # Only accept whitelisted calls for legacy support
        if producer.call.startswith("event::filter::Flag"):
            self.parse_Flag_from_call(
                producer=producer, parent=parent, scope=scope, align=align
            )
        else:
            raise NotImplementedError(
                f"The call {producer.call} does not have legacy support."
            )

    def parse_Flag_from_call(self, producer, parent, scope, align=""):
        """Extract vector configurations from a legacy producer's call string using regex.

        Args:
            producer (object): The legacy flag producer instance.
            parent (str): The ID of the parent node.
            scope (str): The active scope.
            align (str, optional): Debug print alignment spacing. Defaults to "".

        Raises:
            ValueError: If input vector configs cannot be parsed or matched.
            NotImplementedError: If list outputs are used (unsupported in legacy).
        """
        group_id = f"{producer.name}_v"
        call = producer.call
        # Pattern consists of: event::filter::Flag(df, "Flag_name", "Flag_name")
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
        # Determine the index of the input vector config from the call
        vec_input_index = producer.vec_configs.index(input_vec_config)
        call = producer.call
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
        # Loop over all vector configurations
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

            # outputs are not supported for this legacy producer
            if isinstance(producer.output, list):
                raise NotImplementedError(
                    "List outputs for legacy parsed calls are not supported."
                )

    def parse_ExtendedVectorProducer(self, producer, parent, scope, align=""):
        """Parse an ExtendedVectorProducer class instance into graph nodes, inputs, and outputs.

        Args:
            producer (object): The ExtendedVectorProducer instance.
            parent (str): The ID of the parent node.
            scope (str): The active scope.
            align (str, optional): Debug print alignment spacing. Defaults to "".
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

        # ExtendedVectorProducer is better defined and doesn't require regex
        call = producer.call
        vec_config = self.config.config_parameters[scope][producer.vec_config]
        # Loop over all vector configurations
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
        """Parse a ProducerGroup class into graph nodes and recursively route its producers.

        Args:
            producer (object): The ProducerGroup instance.
            parent (str): The ID of the parent node.
            scope (str): The active scope.
            align (str, optional): Debug print alignment spacing. Defaults to "".
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
        # Add all producers in group recursively
        for p in producer.producers[scope]:
            self.parse_Producer_routing(
                producer=p, parent=group_id, scope=scope, align=align + "    "
            )

    def parse_Filter(self, producer, parent, scope, align=""):
        """Parse a Filter class into graph nodes and recursively route its components.

        Args:
            producer (object): The Filter instance.
            parent (str): The ID of the parent node.
            scope (str): The active scope.
            align (str, optional): Debug print alignment spacing. Defaults to "".
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

        # Add all filters in group recursively
        for p in producer.producers[scope]:
            self.parse_Producer_routing(
                producer=p, parent=group_id, scope=scope, align=align + "    "
            )

    def parse_BaseFilter(self, producer, parent, scope, align=""):
        """Parse a legacy BaseFilter class into graph nodes and inputs.

        Args:
            producer (object): The legacy BaseFilter instance.
            parent (str): The ID of the parent node.
            scope (str): The active scope.
            align (str, optional): Debug print alignment spacing. Defaults to "".
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
        # BaseFilter don't have outputs
        for n in set(producer.input[scope]):
            self.add_input(n.name, producer.name, scope)
            if self.debugging:
                print(align + "    " + f"Adding Input: {n.name}")

    def parse_Producer(self, producer, parent, scope, align=""):
        """Parse a generic Producer class into graph nodes, inputs, and outputs.

        Args:
            producer (object): The Producer instance.
            parent (str): The ID of the parent node.
            scope (str): The active scope.
            align (str, optional): Debug print alignment spacing. Defaults to "".
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
        """Designate a node as the source of an Ntuple output.

        Args:
            scope (str): The active scope.
            req_out (str): The name of the requested output quantity.

        Raises:
            ValueError: If the number of producers is invalid, the source is missing,
                        or the output is not provided by NanoAOD/Ntuple.
        """
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
            if not (req_out in self.external_inputs):
                raise ValueError(
                    f"Requested output {req_out} is not provided by NanoAOD/Ntuple."
                )

    def assemble_connections(self):
        """Resolve mappings to construct the actual connecting edges in the DAG.

        Determines the precise source of each input requirement, whether it
        originates from another Producer, NanoAOD, or Ntuple.

        Raises:
            ValueError: If an input is entirely missing from both external sources and internal producers.
        """
        for scope in self.scopes:
            # Assemble connections by iterating through all nodes
            # {target_node: [required_inputs]} is matched to {input: [source_nodes]}
            # Connection is of shape {source: {required_inputs: [targets]}}
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
                    elif req_input in self.external_inputs:
                        if self.is_friend_config:
                            # CROWN friend production may only read quantities from  CROWN Ntuples
                            self.nodes[target_node].setdefault("is_in", set()).add("Ntuple")
                            self.node_family_register[target_node]["file_in"][
                                "Ntuple"
                            ].append(req_input)
                        else:
                            # CROWN Ntuple production may only read quantities from NanoAOD
                            self.nodes[target_node].setdefault("is_in", set()).add(
                                "NanoAOD"
                            )
                            self.node_family_register[target_node]["file_in"][
                                "NanoAOD"
                            ].append(req_input)
                    else:
                        raise ValueError(
                            f"Input {req_input} is missing from NanoAOD/Ntuple and producers."
                        )

                # Create connections grouped by source node
                for source_node, input_names in compose.items():
                    for input_name in input_names:
                        if self.debugging:
                            print(
                                f"Adding Connection {input_name} from {source_node} to {target_node}"
                            )
                        self.add_connection(
                            source=source_node, target=target_node, name=input_name
                        )
        # Convert is_in value to string for cytoscape selector
        for node in self.nodes.values():
            if "is_in" in node:
                node["is_in"] = ",".join(sorted(list(node["is_in"])))

    def bundle_edges(self):
        """Group connections based on common target locations to reduce clutter.

        Creates structural proxy nodes and dynamic edge weighting based on the
        number of overlapping connections.
        """
        # Apply bundling for each quantity separately
        for source, dt in self.connections.items():
            # Inner and outer loop as each source may have multiple output quantities
            # And each quantity may be used by multiple target nodes
            for name, targets in dt.items():
                if len(targets) == 1:
                    # Simplified case: only one target
                    # One proxy node near source for edge label
                    target = targets[0]
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
                    # More than one target: Resolved through relative path bundling
                    # Determine location relative to root
                    source_history = self.get_ancestors(source)
                    # Convert to relative path to source node
                    relative_targets = [
                        self.get_relative(source_history, self.get_ancestors(target))
                        for target in targets
                    ]
                    # bundle based on relative location
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
        """Traverse upwards to return a recursive list of a given node's parents.

        Args:
            node_name (str): The ID of the target node.

        Returns:
            list: A list of parent IDs sequentially ending with the requested node.
        """
        parent = self.nodes[node_name].get("parent")
        if not parent:
            return ["root", node_name]
        else:
            return self.get_ancestors(parent) + [node_name]

    def get_relative(self, source_hist, target_hist):
        """Calculate the relative path from a source tree path to a target tree path.

        Args:
            source_hist (list): The list of ancestors for the source node.
            target_hist (list): The list of ancestors for the target node.

        Returns:
            list: The calculated relative path traversing upwards to the common ancestor and back down.
        """
        min_len = min(len(source_hist), len(target_hist))
        # Find first non-matching element
        for i in range(min_len):
            if source_hist[i] != target_hist[i]:
                # Return shortest path
                return (
                    list(reversed(source_hist[i:]))
                    + [source_hist[i - 1]]
                    + target_hist[i:]
                )

    def construct_proxies(self, name, data, valid_idx, tot_idx, rank, source):
        """Recursively build proxy nodes and edges to visually structure bundled connections.

        Args:
            name (str): The name identifier of the bundled connection.
            data (list): Path history data.
            valid_idx (list): Currently valid indices within the data structure.
            tot_idx (list): Total valid indices at the root operation.
            rank (int): The current depth/rank in the path resolution.
            source (str): The source node to tie the current proxy branch to.
        """
        # Get slice by current rank in path
        data_slice = [data[i][rank] for i in tot_idx if i in valid_idx]
        if len(set(data_slice)) <= 1:
            # Go to next rank if no edge splits off
            self.construct_proxies(name, data, valid_idx, tot_idx, rank + 1, source)
        else:
            # Add a proxy node at the previous (matching) rank
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
            # Run bundling algorithm for each branch recursively
            # Edges now start at proxy
            for val in set(data_slice):
                # Calculate number of connections in the branch
                part_idx = [valid_idx[i] for i, d in enumerate(data_slice) if d == val]
                if len(part_idx) <= 1:
                    # Add leaf edge if only 1 connection remains
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
                    # Go deeper if more than 1 connection remains
                    self.construct_proxies(
                        name, data, part_idx, tot_idx, rank + 1, proxy_node_name
                    )

    def extract_configs(self, call, scope, vector_configs=None, ignore=None):
        """Parse out explicit configuration string replacements from a call parameter.

        Args:
            call (str): The raw call string containing formatting templates.
            scope (str): The active scope context.
            vector_configs (dict, optional): Specific vector configurations to prefer. Defaults to None.
            ignore (set/list, optional): Strings/keys to ignore when extracting. Defaults to None.

        Returns:
            dict: The mapped configurations dict resolving template markers to true values.

        Raises:
            ValueError: If an unknown configuration parameter is encountered.
        """
        # Ignore parameters that are not config parameters
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

        # Get config parameter values from vector or general configs
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
        """Register a node internally as the generator of a specific output.

        Args:
            output_name (str): The name of the output being produced.
            output_node (str): The ID of the node producing it.
            scope (str): The active scope context.

        Raises:
            ValueError: If the node is already registered for this specific output.
        """
        node_id = f"{scope}_{output_node}"
        if node_id in self.outputs[scope][output_name]:
            raise ValueError(f"Node {node_id} already exists in output nodes.")
        self.outputs[scope][output_name].append(node_id)

    def add_input(self, input_name, input_node, scope):
        """Register a node internally as requiring a specific input dependency.

        Args:
            input_name (str): The name of the required input.
            input_node (str): The ID of the node needing it.
            scope (str): The active scope context.

        Raises:
            ValueError: If the input is already registered for this node.
        """
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
        """Create a standardized node object and append it to the graph's internal list.

        Also manages appending metadata to the node and edge family registers.

        Args:
            id_name (str): The base string ID of the node.
            name (str, optional): Visual label for the node. Defaults to id_name.
            scope (str, optional): The active scope context. Defaults to None.
            parent (str, optional): The ID of the parent node. Defaults to None.
            node_type (str, optional): Specific class or type of the node. Defaults to None.
            is_filter (bool, optional): Indicates if the node acts as a filter. Defaults to False.
            family (str, optional): Ties the node to a grouping/proxy family. Defaults to None.
            node_call_data (dict, optional): Stores string call and extraction data. Defaults to None.

        Raises:
            ValueError: If conflicting args (node_type and is_filter) are supplied,
                        or if a full ID already exists within the standard family register.
        """
        if node_type and is_filter:
            raise ValueError("Cannot specify both node_type and is_filter")
        if is_filter:
            node_type = "filter"

        if not name:
            name = id_name
        # Add scope name to node id to distinguish between nodes of the same name
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
        """Log a relational connection locally between a source and a target.

        Args:
            source (str): ID of the source node providing data.
            target (str): ID of the target node consuming data.
            name (str): Connection quantity label.
        """
        self.connections[source][name].append(target)

    def add_edge(self, source, target, edge_id, name, weight, edge_type, family):
        """Create a directional edge drawing parameters for JSON structure export.

        Args:
            source (str): The ID of the node the edge originates from.
            target (str): The ID of the node the edge points towards.
            edge_id (str): The unique ID constructed for this specific edge string.
            name (str): Label describing the edge payload.
            weight (int/float): The raw volume/scale of the edge, formatted later as sqrt(weight).
            edge_type (str): Structural descriptor indicating behavior (e.g., leaf, branch, twig).
            family (str): Tie-in key for grouping in the edge family register.
        """
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

    def save_graph(self, name):
        """Compile the DAG data structures, inject metadata, and export to a JSON file.

        Also triggers an automatic update to the overarching DAG file tracker.

        Args:
            name (str): The base filename used to construct the final exported JSON path.
        """
        path = f"{name}_{self.config.era}_{self.config.sample}_{self.active_scope}.json"
        if self.DAG_dir:
            path = os.path.join(self.DAG_dir, path)
            os.makedirs(self.DAG_dir, exist_ok=True)

        # Compile DAG data with metadata
        meta_data = {
            "nodeFamilyRegister": self.node_family_register,
            "edgeFamilyRegister": self.edge_family_register,
        }
        full_dict = {"elementData": self.DAG_data, "metaData": meta_data}

        # Write DAG data to json
        with open(path, "w") as f:
            json.dump(full_dict, f, indent=4)

        # Update master DAG file list
        self.update_DAG_file_list(
            os.path.join(self.DAG_dir, "DAG_files.json"),
            self.config.era,
            self.config.sample,
            self.active_scope,
        )

    def update_DAG_file_list(self, config_path, new_era, new_sample, new_scope):
        """Maintain and append to the master JSON manifest tracking generated DAG elements.

        Prevents duplicates while verifying the registration of distinct eras, samples, and scopes.

        Args:
            config_path (str): Filepath pointing to the master 'DAG_files.json'.
            new_era (str/int): The era identifier to check/add.
            new_sample (str): The sample identifier to check/add.
            new_scope (str): The scope string to check/add.
        """
        if os.path.exists(config_path):
            with open(config_path, "r") as f:
                try:
                    data = json.load(f)
                except json.JSONDecodeError:
                    data = {"era": [], "sample": [], "scope": []}
        else:
            data = {"era": [], "sample": [], "scope": []}

        # Update master DAG file data
        data["era"] = sorted(list(set(data["era"] + [str(new_era)])))
        data["sample"] = sorted(list(set(data["sample"] + [str(new_sample)])))
        data["scope"] = sorted(list(set(data["scope"] + [str(new_scope)])))

        # Write updated master DAG file data back
        with open(config_path, "w") as f:
            json.dump(data, f, indent=4)

        if self.debugging:
            print(f"Updated {config_path} with {new_era}/{new_sample}/{new_scope}")
