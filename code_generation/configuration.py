from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
from typing import Any, Dict, List, Set, Union

from code_generation.exceptions import (
    ChannelConfigurationError,
    EraConfigurationError,
    InvalidOutputError,
    SampleConfigurationError,
)
from code_generation.modifiers import EraModifier, SampleModifier
from code_generation.optimizer import ProducerOrdering
from code_generation.producer import (
    CollectProducersOutput,
    TProducerInput,
    TProducerStore,
)
from code_generation.quantity import NanoAODQuantity, QuantitiesInput, QuantitiesStore
from code_generation.rules import ProducerRule, RemoveProducer
from code_generation.systematics import SystematicShift, SystematicShiftByQuantity

log = logging.getLogger(__name__)
# type aliases
ScopeSet = Set[str]
ScopeList = Union[str, List[str]]
EraList = Union[str, List[str]]
ChannelList = Union[str, List[str]]
ShiftList = Union[str, List[str]]
SamplesList = Union[str, List[str]]
TConfiguration = Dict[
    str,
    Union[
        List[Union[str, int, float, bool, Dict[Any, Any]]],
        str,
        int,
        float,
        bool,
        EraModifier,
        SampleModifier,
        Dict[Any, Any],
    ],
]


class Configuration(object):
    """
    Configuration class for for the CROWN configuration. This class
    holds all parts of the configuration, from the sample, era, channel,
    and systematics, to the output. All modifications to the configuration should be done through this class.
    """

    def __init__(
        self,
        era: str,
        sample: str,
        channels: ChannelList,
        shifts: ShiftList,
        available_sample_types: SamplesList,
        available_eras: EraList,
        available_channels: ChannelList,
    ):
        """

        Initiation of a new configuration.

        Args:
            era: The era of the sample.
            sample: The sample type of the sample.
            channels: The channels to be used in the configuration. Each channel is considered a scope in CROWN.
            By default, the global scope will always be added and run first, all other scopes will be run in
            parallel as children of the global scope.
            shifts: The systematics to be used in the configuration.
            available_sample_types: The available sample types.
            available_eras: The available eras.
            available_channels: The available channels.

        """
        self.era = era
        self.sample = sample
        self.channels = set(channels)
        self.selected_shifts = set(shifts)
        self.available_sample_types = set(available_sample_types)
        self.available_eras = set(available_eras)
        self.available_channels = set(available_channels)
        self.available_outputs: QuantitiesStore = {}
        self.global_scope = "global"

        self.producers: TProducerStore = {}
        self.scopes: List[str] = []
        self.outputs: QuantitiesStore = {}
        self.shifts: Dict[str, Dict[str, TConfiguration]] = {}
        self.rules: Set[ProducerRule] = set()
        self.config_parameters: Dict[str, TConfiguration] = {}

        self.setup_defaults()

    def validate_channels(self) -> None:
        """
        Function to validate the selected channels. If the channel is not available, an error is raised.

        Args:
            None

        Returns:
            None
        """
        missing_channels = self.channels - self.available_channels
        if len(missing_channels) > 0:
            raise ChannelConfigurationError(missing_channels, self.available_channels)

    def validate_sample_type(self) -> None:
        """
        Function to validate the selected sample type. If the sample type is not available, an error is raised.

        Args:
            None

        Returns:
            None
        """
        if self.sample not in self.available_sample_types:
            raise SampleConfigurationError(self.sample, self.available_sample_types)

    def validate_era(self) -> None:
        """
        Function to validate the selected era. If the era is not available, an error is raised.

        Args:
            None

        Returns:
            None
        """
        if self.era not in self.available_eras:
            raise EraConfigurationError(self.era, self.available_eras)

    def set_sample_parameters(self) -> None:
        """
        Helper function to add sample type variables to the configuration.
        The variables look like "is_${sampletype}" and can be used in all producer
        calls to check the type of the sample.

        Args:
            None

        Returns:
            None
        """
        sample_parameters: Dict[str, bool] = {}
        for sampletype in self.available_sample_types:
            if self.sample == sampletype:
                sample_parameters["is_{}".format(sampletype)] = True
            else:
                sample_parameters["is_{}".format(sampletype)] = False
        for scope in self.scopes:
            self.config_parameters[scope].update(sample_parameters)

    def setup_defaults(self) -> None:
        """
        Function used to add some defaults to the configuration. This function is called by the __init__ function.
        For all configured channels, the nessessay variables are added to the configuration.
        The validation of the initial settings is also done here.

        Args:
            None

        Returns:
            None
        """
        self.validate_channels()
        self.validate_sample_type()
        self.validate_era()
        self.scopes = [self.global_scope]
        for channel in self.available_channels:
            self.scopes.append(channel)
        for scope in self.scopes:
            self.producers[scope] = []
            self.outputs[scope] = set()
            self.shifts[scope] = {}
            self.available_outputs[scope] = set()
            self.config_parameters[scope] = {}
        self.set_sample_parameters()

    def add_config_parameters(
        self, scopes: ScopeList, parameters: TConfiguration
    ) -> None:
        """
        Function to add new config parameters to the configuration. Modifiers are used to
        modify the configuration, are directly resolved using the resolve_modifiers function.

        Args:
            scopes: The scopes to which the parameters should be added. This can be a list of scopes or a single scope.
            parameters: The parameters to be added. This must be a dictionary of parameters. If multiple scopes are given,
            the parameters are added to all scopes.

        Returns:
            None
        """
        if not isinstance(scopes, list):
            scopes = [scopes]
        for scope in scopes:
            self.config_parameters[scope].update(self.resolve_modifiers(parameters))

    def add_producers(self, scopes: ScopeList, producers: TProducerInput) -> None:
        """
        Function used to add producers to the configuration. Internally, a set of all
        available outputs is updated, which is later used to check if all required ouputs are available.

        Args:
            scopes: The scopes to which the producers should be added. This can be a list of scopes or a single scope.
            producers: The producers to be added. This must be a list of producers. If multiple scopes are given,
            the producers are added to all scopes.

        Returns:
            None
        """
        if not isinstance(scopes, list):
            scopes = [scopes]
        if not isinstance(producers, list):
            producers = [producers]
        for scope in scopes:
            self.producers[scope].extend(producers)
            self.available_outputs[scope].update(
                CollectProducersOutput(producers, scope)
            )

    def add_outputs(self, scopes: ScopeList, output: QuantitiesInput) -> None:
        """
        Function used to add outputs to the configuration.

        Args:
            scopes: The scopes to which the outputs should be added. This can be a list of scopes or a single scope.
            outputs: The outputs to be added. This must be a list of outputs. If multiple scopes are given,
            the outputs are added to all scopes.

        Returns:
            None
        """
        if not isinstance(scopes, list):
            scopes = [scopes]
        if not isinstance(output, list):
            output = [output]
        for scope in scopes:
            self.outputs[scope].update(output)

    def add_shift(
        self,
        shift: Union[SystematicShift, SystematicShiftByQuantity],
        samples: Union[str, List[str], None] = None,
    ) -> None:
        """
        Function used to add a systematics shift to the configuration. During this step, the shift is validated and applied.

        Args:
            shift: The shift to be added. This must be a SystematicShift object.
            samples: The samples to which the shift should be applied. This can be a list of samples or a single sample.
                If ths option is not set, the shift is applied, regardless of the sample type.

        Returns:
            None
        """
        if not isinstance(shift, SystematicShift):
            raise TypeError("shift must be of type SystematicShift")
        if isinstance(samples, str):
            samples = [samples]
        if samples is None or self.sample in samples:
            scopes_to_shift: ScopeList = [
                scope for scope in shift.get_scopes() if scope in self.scopes
            ]
            for scope in scopes_to_shift:
                shift.apply(scope)
                config_change = shift.get_shift_config(scope)
                shiftname = shift.shiftname
                self.shifts[scope][shiftname] = config_change

    def add_modification_rule(self, scopes: ScopeList, rule: ProducerRule) -> None:
        """
        Function used to add a rule to the configuration.

        Args:
            scopes: The scopes to which the rule should be added. This can be a list of scopes or a single scope.
            rule: The rule to be added. This must be a ProducerRule object.

        Returns:
            None

        """
        if not isinstance(rule, ProducerRule):
            raise TypeError("Rule must be of type ProducerRule")
        if not isinstance(scopes, list):
            scopes = [scopes]
        rule.set_scopes(scopes)
        rule.set_global_scope(self.global_scope)
        self.rules.add(rule)

    def resolve_modifiers(self, configuration_dict: Dict[Any, Any]) -> TConfiguration:
        """
        Function used to resolve mofifiers used in the configuration. This function is called by the add_config_parameters function.

        Args:
            configuration_dict: The configuration dictionary to be resolved.

        Returns:
            dict: The resolved configuration dictionary.
        """
        resolved_dict: TConfiguration = {}
        for key in configuration_dict:
            resolved_value = None
            log.debug("Resolved: {}".format(resolved_dict))
            log.debug("Parameter: {}".format(key))
            log.debug("Paramters: {}".format(configuration_dict))
            if isinstance(configuration_dict[key], dict):
                resolved_dict.update(self.resolve_modifiers(configuration_dict[key]))
            elif isinstance(configuration_dict[key], SampleModifier):
                log.debug("Applying Samplemodifier")
                resolved_value = configuration_dict[key].apply(self.sample)
            elif isinstance(configuration_dict[key], EraModifier):
                log.debug("Applying Eramodifier")
                resolved_value = configuration_dict[key].apply(self.era)
            else:
                resolved_value = configuration_dict[key]
            resolved_dict.update({key: resolved_value})
        return resolved_dict

    def remove_empty_scopes(self) -> None:
        """
        Function used to remove empty scopes from the configuration. This function is called by the optimize function,
        which should be called after all configuration changes have been made.

        Args:
            None

        Returns:
            None
        """
        # we have to use a seperate list, because we cannot modify the list while iterating over it without breaking stuff
        scopes_to_test = [scope for scope in self.scopes]
        for scope in scopes_to_test:
            if (len(self.producers[scope]) == 0) or (
                scope not in self.channels and scope is not self.global_scope
            ):
                log.warning("Removing unrequested / empty scope {}".format(scope))
                self.scopes.remove(scope)
                del self.producers[scope]
                del self.outputs[scope]
                del self.shifts[scope]
                del self.config_parameters[scope]
                del self.available_outputs[scope]

    def apply_rules(self) -> None:
        """
        Function used to apply all rules to the configuration. This function is called by the optimize function.

        Args:
            None

        Returns:
            None
        """
        for rule in self.rules:
            rule.apply(self.sample, self.producers, self.outputs)
            # also update the set of available outputs
            for scope in rule.affected_scopes():
                if isinstance(rule, RemoveProducer):
                    self.available_outputs[scope] - CollectProducersOutput(
                        rule.affected_producers(), scope
                    )
                else:
                    self.available_outputs[scope].update(
                        CollectProducersOutput(rule.affected_producers(), scope)
                    )

    def optimize(self) -> None:
        """
        Function used to optimize the configuration. Optimizaion steps are:
        - Remove empty scopes
        - Apply rules
        - Optimizing producer ordering (this does not change the configuration, but only the order of producers)

        Args:
            None

        Returns:
            None
        """
        self.apply_rules()
        self.remove_empty_scopes()
        for scope in self.scopes:
            log.debug("Optimizing Producer Ordering in scope {}".format(scope))
            ordering = ProducerOrdering(
                global_producers=list(self.producers[self.global_scope]),
                scope=scope,
                producer_ordering=list(self.producers[scope]),
            )
            ordering.Optimize()
            self.producers[scope] = ordering.optimized_ordering

    def validate_outputs(self) -> None:
        """
        Function used to validate the defined outputs. If an output is requested in the configuration,
        but is not available, since no producer will be able to produce it, an error is raised.

        Args:
            None

        Returns:
            None
        """
        for scope in [scope for scope in self.scopes if scope != self.global_scope]:
            required_outputs = set(
                output
                for output in self.outputs[scope] | self.outputs[self.global_scope]
                if not isinstance(output, NanoAODQuantity)
            )
            # merge the two sets of outputs
            provided_outputs = (
                self.available_outputs[scope]
                | self.available_outputs[self.global_scope]
            )
            missing_outputs = required_outputs - provided_outputs
            if len(missing_outputs) > 0:
                raise InvalidOutputError(scope, missing_outputs)

    def validate(self) -> None:
        """
        Function used to validate the configuration. During the validation, the following steps are performed:
        - Validate the outputs

        Args:
            None

        Returns:
            None
        """
        self.validate_outputs()

    def report(self) -> None:
        """
        Function used to log a summary of the configuration in form of a report.

        Args:
            None

        Returns:
            None
        """
        running_scopes = list(self.channels) + [self.global_scope]
        total_producers = sum([len(self.producers[scope]) for scope in running_scopes])
        total_quantities = sum([len(self.outputs[scope]) for scope in running_scopes])
        total_shifts = sum([len(self.shifts[scope]) for scope in running_scopes])
        log.info("------------------------------------")
        log.info("Configuration Report")
        log.info("------------------------------------")
        log.info("  Sample: {}".format(self.sample))
        log.info("  Era: {}".format(self.era))
        log.info("  Channels: {}".format(self.channels))

        log.info("  Total number of producers: {}".format(total_producers))
        for scope in running_scopes:
            log.info("       {}: {}".format(scope, len(self.producers[scope])))
        log.info("  Total number of quantities: {}".format(total_quantities))
        for scope in running_scopes:
            log.info("       {}: {}".format(scope, len(self.outputs[scope])))
        log.info("  Total number of shifts: {}".format(total_shifts))
        for scope in running_scopes:
            log.info("       {}: {}".format(scope, len(self.shifts[scope])))
        log.info("------------------------------------")

    def __str__(self) -> str:
        """
        String representation of the configuration.
        """
        returnstr = "Configuration:"
        returnstr += "  Era: {}".format(self.era)
        returnstr += "  Sample: {}".format(self.sample)
        returnstr += "  Channels: {}".format(self.channels)
        returnstr += "  Shifts: {}".format(self.shifts)
        returnstr += "  Rules:  {}".format(self.rules)
        returnstr += "  Outputs:"
        for scope in self.scopes:
            returnstr += "    {}: {}".format(scope, self.outputs[scope])
        returnstr += "  Producers:"
        for scope in self.scopes:
            returnstr += "    {}: {}".format(scope, self.producers[scope])
        returnstr += "  Config Parameters:"
        for scope in self.scopes:
            returnstr += "    {}: {}".format(scope, self.config_parameters[scope])
        return returnstr

    def dump_dict(self) -> Dict[str, Dict[Any, Any]]:
        """
        Function used to dump the configuration as a dict, which is used to generate the c++ file.

        Args:
            None

        Returns:
            dict: The configuration as a dict
        """
        returndict: Dict[str, Dict[Any, Any]] = {}
        returndict[""] = {}
        returndict["output"] = {}
        returndict["producers"] = {}
        for scope in self.scopes + [self.global_scope]:
            returndict[""][scope] = self.config_parameters[scope]
            returndict["producers"][scope] = self.producers[scope]
            if scope is not self.global_scope:
                log.debug(
                    "Final set of outputs : {}".format(
                        sorted(list(self.outputs[scope]))
                    )
                )
                returndict["output"][scope] = sorted(list(self.outputs[scope]))
            # add systematic shifts
            for shift in self.shifts[scope]:
                log.warning("Adding shift {} in scope {}".format(shift, scope))
                log.debug("  {}".format(self.shifts[scope][shift]))
                try:
                    returndict[shift][scope] = (
                        self.config_parameters[scope] | self.shifts[scope][shift]
                    )
                except KeyError:
                    returndict[shift] = {}
                    returndict[shift][scope] = (
                        self.config_parameters[scope] | self.shifts[scope][shift]
                    )
        return returndict
