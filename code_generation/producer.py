from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
from typing import Any, Dict, List, Set, Union

from code_generation.exceptions import (
    InvalidProducerConfigurationError,
    ConfigurationError,
)

import code_generation.quantity as q

log = logging.getLogger(__name__)


class SafeDict(Dict[Any, Any]):
    def __missing__(self, key: Any) -> Any:
        return "{" + key + "}"


class Producer:
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        output: Union[List[q.Quantity], None],
        scopes: List[str],
    ):
        """Docstrings to be added ...."""
        log.debug("Setting up a new producer {}".format(name))

        # sanity checks
        if not isinstance(input, list) and not isinstance(input, dict):
            log.error(
                "Exception (%s): Argument 'input' must be a list or a dict!" % name
            )
            raise Exception
        if not isinstance(output, list) and output is not None:
            log.error(
                "Exception (%s): Argument 'output' must be a list or None!" % name
            )
            raise Exception
        self.name: str = name
        self.call: str = call
        self.output: Union[List[q.Quantity], None] = output
        self.scopes = scopes
        # if input not given as dict and therfore not scope specific transform into dict with all scopes
        if not isinstance(input, dict):
            inputdict = {}
            for scope in self.scopes:
                inputdict[scope] = input.copy() if isinstance(input, list) else input
        else:
            inputdict = input
        self.input: Dict[str, List[q.Quantity]] = inputdict
        # keep track of variable dependencies
        if self.output is not None:
            for scope in self.scopes:
                for input_quantity in self.input[scope]:
                    for output_quantity in self.output:
                        input_quantity.adopt(output_quantity, scope)
        log.debug("-----------------------------------------")
        log.debug("| Producer: {}".format(self.name))
        log.debug("| Call: {}".format(self.call))
        for scope in self.scopes:
            if self.input[scope] is None:
                log.debug("| Inputs ({}): None".format(scope))
            else:
                log.debug(
                    "| Inputs ({}): {}".format(
                        scope, [input.name for input in self.input[scope]]
                    )
                )
        if self.output is None:
            log.debug("| Output: None")
        else:
            log.debug("| Outputs: {}".format([output.name for output in self.output]))
        log.debug("| scopes: {}".format(self.scopes))
        log.debug("-----------------------------------------")

    def __str__(self) -> str:
        return "Producer: {}".format(self.name)

    def __repr__(self) -> str:
        return "Producer: {}".format(self.name)

    # Check if a output_quantity is already used as an output by
    # another producer within the same scope.
    # If this occurs, a Exception is thrown, since this is not possible with dataframes
    def reserve_output(self, scope: str) -> None:
        if self.output is not None:
            for output_quantity in self.output:
                output_quantity.reserve_scope(scope)

    def shift(self, name: str, scope: str = "global") -> None:
        if scope not in self.scopes:
            log.error(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        if self.output is None:
            log.error(
                "Exception (%s): output None cannot be shifted ! How did you end up here ?"
                % name
            )
            raise Exception
        for entry in self.output:
            entry.shift(name, scope)

    def ignore_shift(self, name: str, scope: str = "global") -> None:
        if scope not in self.scopes:
            log.error(
                "Trying to add shift %s to producer %s in scope %s, but producer does not exist in given scope!"
                % (name, self.name, scope)
            )
            raise Exception
        if self.output is None:
            log.error(
                "Exception (%s): output None cannot be shifted ! How did you end up here ?"
                % name
            )
            raise Exception
        for entry in self.output:
            entry.ignore_shift(name, scope)

    def writecall(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str, shift: str = ""
    ) -> str:
        if self.output is None:
            config[shift][scope]["output"] = ""
            config[shift][scope]["output_vec"] = ""
        else:
            log.debug("Writing call for {}".format(self.name))
            log.debug("output: {}".format(self.output))
            log.debug("name: {}".format(self.name))
            log.debug("shift: {}".format(shift))
            log.debug("Scopes: {}".format(config[shift].keys()))
            config[shift][scope]["output"] = (
                '"' + '","'.join([x.get_leaf(shift, scope) for x in self.output]) + '"'
            )
            config[shift][scope]["output_vec"] = (
                '{"'
                + '","'.join([x.get_leaf(shift, scope) for x in self.output])
                + '"}'
            )
        config[shift][scope]["input"] = (
            '"'
            + '", "'.join([x.get_leaf(shift, scope) for x in self.input[scope]])
            + '"'
        )
        config[shift][scope]["input_vec"] = (
            '{"'
            + '","'.join([x.get_leaf(shift, scope) for x in self.input[scope]])
            + '"}'
        )
        config[shift][scope]["df"] = "{df}"
        config[shift][scope]["vec_open"] = "{vec_open}"
        config[shift][scope]["vec_close"] = "{vec_close}"
        # if a bool is used in the python configuration, convert it to a c++ bool value
        # True -> true, False -> false
        for para in config[shift][scope]:
            if isinstance(config[shift][scope][para], bool):
                if config[shift][scope][para]:
                    log.debug("Found a boolean True ! - converting to C++ syntax")
                    config[shift][scope][para] = "true"
                else:
                    log.debug("Found a boolean False ! - converting to C++ syntax")
                    config[shift][scope][para] = "false"
        try:
            return self.call.format(
                **config[shift][scope]
            )  # use format (not format_map here) such that missing config entries cause an error
        except KeyError as e:
            log.error(
                "Error in {} Producer, key {} is not found in configuration".format(
                    self.name, e
                )
            )
            log.error(config[shift][scope])
            log.error("Call: {}".format(self.call))
            raise Exception

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str
    ) -> List[str]:
        if scope not in self.scopes:
            log.error(
                "Exception ({}): Tried to use producer in scope {}, which the producer is not forseen for!".format(
                    self.name, scope
                )
            )
            raise Exception
        calls = [self.writecall(config, scope)]
        if self.output is not None:
            list_of_shifts = self.output[0].get_shifts(
                scope
            )  # all entries must have same shifts
            for shift in list_of_shifts:
                calls.append(self.writecall(config, scope, shift))
        return calls

    def get_inputs(self, scope: str) -> List[q.Quantity]:
        if scope not in self.scopes:
            log.error(
                "Exception ({}): Tried to get producer inputs in scope {}, which the producer is not forseen for!".format(
                    self.name, scope
                )
            )
            raise Exception
        return self.input[scope]

    def get_outputs(self, scope: str) -> List[Union[q.QuantityGroup, q.Quantity]]:
        if scope not in self.scopes:
            log.error(
                "Exception ({}): Tried to get producer outputs in scope {}, which the producer is not forseen for!".format(
                    self.name, scope
                )
            )
            raise Exception
        if self.output is None:
            return []
        else:
            return self.output


class VectorProducer(Producer):
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        output: Union[List[q.Quantity], None],
        scopes: List[str],
        vec_configs: List[str],
    ):
        self.name = name
        super().__init__(name, call, input, output, scopes)
        self.vec_configs = vec_configs

    def __str__(self) -> str:
        return "VectorProducer: {}".format(self.name)

    def __repr__(self) -> str:
        return "VectorProducer: {}".format(self.name)

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str
    ) -> List[str]:
        basecall = self.call
        calls: List[str] = []
        shifts = [""]
        if self.output is not None:
            shifts.extend(self.output[0].get_shifts(scope))
        for shift in shifts:
            # check that all config lists (and output if applicable) have same length
            n_versions = len(config[shift][scope][self.vec_configs[0]])
            for key in self.vec_configs:
                if n_versions != len(config[shift][scope][key]):
                    log.error(
                        "Following lists in config must have same length: %s, %s"
                        % (self.vec_configs[0], key)
                    )
                    raise Exception
            if self.output is not None and len(self.output) != n_versions:
                log.error(
                    "{} expects either no output or same amount as entries in config lists !".format(
                        self
                    )
                )
                log.error("Number of expected outputs: {}".format(n_versions))
                log.error("List of outputs: {}".format(self.output))
                raise Exception
            for i in range(n_versions):
                helper_dict: Dict[Any, Any] = {}
                for key in self.vec_configs:
                    helper_dict[key] = config[shift][scope][key][i]
                if self.output is not None:
                    helper_dict["output"] = (
                        '"' + self.output[i].get_leaf(shift, scope) + '"'
                    )
                    helper_dict["output_vec"] = (
                        '{"' + self.output[i].get_leaf(shift, scope) + '"}'
                    )
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, scope, shift))
        self.call = basecall
        return calls


class ExtendedVectorProducer(Producer):
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        output: str,
        scope: Union[List[str], str],
        vec_config: str,
    ):
        # we create a Quantity Group, which is updated during the writecalls() step
        self.outputname = output
        self.vec_config = vec_config
        if not isinstance(scope, list):
            scope = [scope]
        if len(scope) != 1:
            log.error("ExtendedVectorProducer can only use one scope per instance !")
            raise Exception
        quantity_group = q.QuantityGroup(name)
        # set the vec config key of the quantity group
        quantity_group.set_vec_config(vec_config)
        super().__init__(name, call, input, [quantity_group], scope)
        if self.output is None:
            raise InvalidProducerConfigurationError(self.name)

    def __str__(self) -> str:
        return "ExtendedVectorProducer: {}".format(self.name)

    def __repr__(self) -> str:
        return "ExtendedVectorProducer: {}".format(self.name)

    @property
    def output_group(self) -> q.QuantityGroup:
        if self.output is None:
            raise Exception("ExtendedVectorProducer has no output!")
        if not isinstance(self.output[0], q.QuantityGroup):
            log.error("ExtendedVectorProducer expects a QuantityGroup as output!")
            raise Exception
        return self.output[0]

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, Any]]], scope: str
    ) -> List[str]:
        n_versions = len(config[""][scope][self.vec_config])
        log.debug("Number of extended producers to be created {}".format(n_versions))
        if self.output is None:
            raise InvalidProducerConfigurationError(self.name)
        if not isinstance(self.output[0], q.QuantityGroup):
            log.error("ExtendedVectorProducer expects a QuantityGroup as output!")
            raise Exception
        for i in range(n_versions):
            self.output[0].add(config[""][scope][self.vec_config][i][self.outputname])
        basecall = self.call
        calls: List[str] = []
        shifts = [""]
        shifts.extend(self.output[0].get_shifts(scope))
        for shift in shifts:
            for i in range(n_versions):
                # the information for the producer is directly read from the configuration
                helper_dict = config[shift][scope][self.vec_config][i]
                helper_dict["output"] = (
                    '"' + self.output[0].quantities[i].get_leaf(shift, scope) + '"'
                )
                helper_dict["output_vec"] = (
                    '{"' + self.output[0].quantities[i].get_leaf(shift, scope) + '"}'
                )
                self.call = basecall.format_map(SafeDict(helper_dict))
                calls.append(self.writecall(config, scope, shift))
        self.call = basecall
        return calls


class BaseFilter(Producer):
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        scopes: List[str],
    ):
        super().__init__(name, call, input, None, scopes)

    def __str__(self) -> str:
        return "BaseFilter: {}".format(self.name)

    def __repr__(self) -> str:
        return "BaseFilter: {}".format(self.name)

    def writecall(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str, shift: str = ""
    ) -> str:
        log.critical("{}: Filters do not support method writecall!".format(self.name))
        raise Exception

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str
    ) -> List[str]:
        inputs: List[str] = []
        for quantity in self.input[scope]:
            inputs.extend(quantity.get_leaves_of_scope(scope))
        config[""][scope]["input"] = '"' + '", "'.join(inputs) + '"'
        config[""][scope]["input_vec"] = '{"' + '","'.join(inputs) + '"}'
        config[""][scope]["df"] = "{df}"
        try:
            return [
                self.call.format(**config[""][scope])
            ]  # use format (not format_map here) such that missing config entries cause an error
        except KeyError as e:
            log.error(
                "Error in {} Basefilter, key {} is not found in configuration".format(
                    self.name, e
                )
            )
            log.error("Call: {}".format(self.call))
            raise Exception


class ProducerGroup:
    PG_count = 1  # counter for internal quantities used by ProducerGroups

    def __init__(
        self,
        name: str,
        call: Union[str, None],
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]], None],
        output: Union[List[q.Quantity], None],
        scopes: List[str],
        subproducers: Union[
            List[Producer | ProducerGroup],
            Dict[str, List[Producer | ProducerGroup]],
        ],
    ):
        self.name = name
        self.call = call
        self.output = output
        self.scopes = scopes
        self.producers: Dict[str, List[Producer | ProducerGroup]] = {}
        # if subproducers are given as dict and therefore scope specific transform into dict with all scopes
        if not isinstance(subproducers, dict):
            log.debug("Converting subproducer list to dictionary")
            for scope in self.scopes:
                self.producers[scope] = [producer for producer in subproducers]
        else:
            self.producers = subproducers
        # do a consistency check for the scopes
        self.check_producer_scopes()
        # if input not given as dict and therefore not scope specific transform into dict with all scopes
        if not isinstance(input, dict):
            inputdict = {}
            for scope in self.scopes:
                inputdict[scope] = input.copy() if isinstance(input, list) else input
            self.input = inputdict
        else:
            self.input = dict(input)
        # If call is provided, this is supposed to consume output of subproducers. Creating these internal products below:
        if self.call is not None:
            log.debug("Constructing {}".format(self.name))
            log.debug(" --> Scopes: {}".format(self.scopes))
            for scope in self.scopes:
                log.debug("Adding for scope: {}".format(scope))
                log.debug(" --> Producers {}".format(self.producers[scope]))
                for subproducer in self.producers[scope]:
                    # skip producers without output
                    log.debug(
                        "Adding {} / {} : output: {}".format(
                            subproducer, scope, subproducer.output
                        )
                    )
                    if subproducer.output is None:
                        continue
                    # if the subproducer does not have an output quantity, we assign an internally tracked quantity
                    if subproducer.output == []:
                        # create quantities that are produced by subproducers and then collected by the final call of the producer group
                        subproducer.output = [
                            q.Quantity(
                                "PG_internal_quantity_%i" % self.__class__.PG_count
                            )
                        ]  # quantities of vector producers will be duplicated later on when config is known
                        log.debug(
                            "Added new internal quantitiy: {}".format(
                                subproducer.output
                            )
                        )
                        self.__class__.PG_count += 1
                        if isinstance(subproducer, ProducerGroup):
                            subproducer.producers[scope][-1].output = subproducer.output
                        for subproducer_scope in subproducer.scopes:
                            for quantity in subproducer.input[subproducer_scope]:
                                quantity.adopt(subproducer.output[0], subproducer_scope)
                    for output_quantity in subproducer.output:
                        if output_quantity not in self.input[scope]:
                            self.input[scope].append(output_quantity)
            # treat own collection function as subproducer
            self.setup_own_producer()
        log.debug("-----------------------------------------")
        log.debug("| ProducerGroup: {}".format(self.name))
        log.debug("| Call: {}".format(self.call))
        for scope in self.scopes:
            if self.input[scope] is None:
                log.debug("| Inputs ({}): None".format(scope))
            else:
                log.debug(
                    "| Inputs ({}): {}".format(
                        scope, [input.name for input in self.input[scope]]
                    )
                )
        if self.output is None:
            log.debug("| Output: None")
        else:
            log.debug("| Outputs: {}".format([output.name for output in self.output]))
        log.debug("| scopes: {}".format(self.scopes))
        for scope in self.scopes:
            log.debug(
                "| Producers ({}): {}".format(
                    scope, [producer.name for producer in self.producers[scope]]
                )
            )

        log.debug("-----------------------------------------")

    def __str__(self) -> str:
        return "ProducerGroup: {}".format(self.name)

    def __repr__(self) -> str:
        return "ProducerGroup: {}".format(self.name)

    def check_producer_scopes(self) -> None:
        for scope in self.scopes:
            if scope not in self.producers.keys():
                raise Exception(
                    "ProducerGroup {}: scope {} not found in subproducer configuration: {}".format(
                        self.name, scope, self.producers
                    )
                )

    def setup_own_producer(self) -> None:
        producer = Producer(self.name, self.call, self.input, self.output, self.scopes)
        for scope in self.scopes:
            self.producers[scope].append(producer)

    # for a producer group, step iteratively
    # through the subproducers and reserve the output there
    def reserve_output(self, scope: str) -> None:
        for subproducer in self.producers[scope]:
            subproducer.reserve_output(scope)

    def shift(self, name: str, scope: str = "global") -> None:
        for producer in self.producers[scope]:
            producer.shift(name, scope)

    def ignore_shift(self, name: str, scope: str = "global") -> None:
        for producer in self.producers[scope]:
            producer.ignore_shift(name, scope)

    def writecall(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str, shift: str = ""
    ) -> str:
        raise NotImplementedError("This function is not supported to a ProducerGroup!")

    def writecalls(
        self, config: Dict[str, Dict[str, Dict[str, str]]], scope: str
    ) -> List[str]:
        calls: List[str] = []
        for producer in self.producers[scope]:
            # duplicate outputs of vector subproducers if they were generated automatically
            if (
                self.call is not None
                and isinstance(producer, VectorProducer)
                and producer.output is not None
            ):
                for i in range(len(config[""][scope][producer.vec_configs[0]]) - 1):
                    producer.output.append(
                        producer.output[0].copy(
                            "PG_internal_quantity_%i" % self.__class__.PG_count
                        )
                    )
                    self.__class__.PG_count += 1
                    if producer.output[-1] not in self.input[scope]:
                        self.input[scope].append(producer.output[-1])
            # retrieve calls of subproducers
            calls.extend(producer.writecalls(config, scope))
        return calls

    def get_inputs(self, scope: str) -> List[q.Quantity]:
        inputs: List[q.Quantity] = []
        log.debug("Getting inputs for {}".format(self.name))
        for subproducer in self.producers[scope]:
            log.debug("  --> {} {}".format(subproducer, subproducer.get_inputs(scope)))
            inputs.extend(subproducer.get_inputs(scope))
        return inputs

    def get_outputs(self, scope: str) -> List[q.Quantity]:
        outputs: List[q.Quantity] = []
        log.debug("Getting outputs for {}".format(self.name))
        for subproducer in self.producers[scope]:
            log.debug("  --> {} {}".format(subproducer, subproducer.get_outputs(scope)))
            outputs.extend(subproducer.get_outputs(scope))
        return outputs


class Filter(ProducerGroup):
    def __init__(
        self,
        name: str,
        call: str,
        input: Union[List[q.Quantity], Dict[str, List[q.Quantity]]],
        scopes: List[str],
        subproducers: Union[
            List[Producer | ProducerGroup],
            Dict[str, List[Producer | ProducerGroup]],
        ],
    ):
        self.__class__.PG_count = ProducerGroup.PG_count
        super().__init__(name, call, input, None, scopes, subproducers)
        ProducerGroup.PG_count = self.__class__.PG_count
        self.call: str = call

    def __str__(self) -> str:
        return "Filter: {}".format(self.name)

    def __repr__(self) -> str:
        return "Filter: {}".format(self.name)

    def setup_own_producer(self) -> None:
        for scope in self.scopes:
            self.producers[scope].append(
                BaseFilter(self.name, self.call, self.input, self.scopes)
            )


def CollectProducersOutput(
    producers: List[Union[ProducerGroup, Producer]], scope: str
) -> Set[q.Quantity]:
    output: Set[q.Quantity] = set()
    for producer in producers:
        if producer.output is not None:
            output |= set(producer.output)
        if isinstance(producer, ProducerGroup):
            try:
                for prod in producer.producers[scope]:
                    output |= CollectProducerOutput(prod, scope)
            except KeyError:
                raise ConfigurationError(
                    "Scope {} not found in subproducer configuration: {}".format(
                        scope, producer.producers
                    )
                )
    return set(output)


def CollectProducerOutput(
    producer: Union[ProducerGroup, Producer], scope: str
) -> Set[q.Quantity]:
    output: Set[q.Quantity] = set()
    if producer.output is not None:
        output |= set(producer.output)
    if isinstance(producer, ProducerGroup):
        try:
            for prod in producer.producers[scope]:
                output |= CollectProducerOutput(prod, scope)
        except KeyError:
            raise ConfigurationError(
                "Scope {} not found in subproducer configuration: {}".format(
                    scope, producer.producers
                )
            )
    return set(output)


TProducerInput = Union[Producer, ProducerGroup, List[Union[Producer, ProducerGroup]]]
TProducerSet = Union[Producer, ProducerGroup, Set[Union[Producer, ProducerGroup]]]
TProducerStore = Dict[str, List[Union[Producer, ProducerGroup]]]
TProducerListStore = Dict[str, TProducerStore]
