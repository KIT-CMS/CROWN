from code_generation.exceptions import (
    ChannelConfigurationError,
    SampleConfigurationError,
    EraConfigurationError,
    InvalidOutputError,
)
from config.utility import CollectProducersOutput
from code_generation.optimizer import ProducerOrdering
from code_generation.quantity import NanoAODQuantity
from code_generation.systematics import SystematicShift
import logging


log = logging.getLogger(__name__)

"""
    Configuration class for for the CROWN configuration. This class holds all parts of the configuration, from the sample, era, channel, and systematics, to the output. All modifications to the configuration should be done through this class.
"""


class Configuration(object):
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

    def __init__(
        self,
        era,
        sample,
        channels,
        shifts,
        available_sample_types,
        available_eras,
        available_channels,
    ):
        self.era = era
        self.sample = sample
        self.channels = set(channels)
        self.shifts = set(shifts)
        self.available_sample_types = set(available_sample_types)
        self.available_eras = set(available_eras)
        self.available_channels = set(available_channels)
        self.available_outputs = {}
        self.global_scope = "global"

        self.producers = {}
        self.scopes = []
        self.outputs = {}
        self.shifts = {}
        self.rules = set()
        self.config_parameters = {}

        self.setup_defaults()

    """
    Function to validate the selected channels. If the channel is not available, an error is raised.

    Args:
        None 
    Returns:
        None  
    """

    def validate_channels(self):
        missing_channels = self.channels - self.available_channels
        if len(missing_channels) > 0:
            raise ChannelConfigurationError(missing_channels)

    """
    Function to validate the selected sample type. If the sample type is not available, an error is raised.

    Args:
        None 
    Returns:
        None  
    """

    def validate_sample_type(self):
        if self.sample not in self.available_sample_types:
            raise SampleConfigurationError(self.sample, self.available_sample_types)

    """
    Function to validate the selected era. If the era is not available, an error is raised.

    Args:
        None 
    Returns:
        None  
    """

    def validate_era(self):
        if self.era not in self.available_eras:
            raise EraConfigurationError(self.era, self.available_eras)

    """
    Helper function to add sample type variables to the configuration. The variables look like 
    "is_${sampletype}" and can be used in all producer calls to check the type of the sample.

    Args:
        None
    Returns:
        None
    """

    def set_sample_parameters(self):
        sample_parameters = {}
        for sampletype in self.available_sample_types:
            if self.sample == sampletype:
                sample_parameters["is_{}".format(sampletype)] = True
            else:
                sample_parameters["is_{}".format(sampletype)] = False
        for scope in self.scopes:
            self.config_parameters[scope].update(sample_parameters)

    """
    Function used to add some defaults to the configuration. This function is called by the __init__ function.
    For all configured channels, the nessessay variables are added to the configuration.
    The validation of the initial settings is also done here.

    Args:  
        None
    Returns:
        None
    """

    def setup_defaults(self):
        self.validate_channels()
        self.validate_sample_type()
        self.validate_era()
        self.scopes = [self.global_scope]
        for channel in self.available_channels:
            self.scopes.append(channel)
        for scope in self.scopes:
            self.producers[scope] = set()
            self.outputs[scope] = set()
            self.shifts[scope] = {}
            self.available_outputs[scope] = set()
            self.config_parameters[scope] = {}
        self.set_sample_parameters()

    """
    Function to add new config parameters to the configuration. 

    Args:
        scopes: The scopes to which the parameters should be added. This can be a list of scopes or a single scope.
        parameters: The parameters to be added. This must be a dictionary of parameters. If multiple scopes are given, 
        the parameters are added to all scopes.
    Returns:
        None
    """

    def add_config_parameters(self, scopes, parameters):
        if not isinstance(scopes, list):
            scopes = [scopes]
        for scope in scopes:
            self.config_parameters[scope].update(self.resolve_modifiers(parameters))

    """
    Function used to add producers to the configuration. Internally, a set of all available outputs is updated, which is later used to check if all required ouputs are available.

    Args:
        scopes: The scopes to which the producers should be added. This can be a list of scopes or a single scope.
        producers: The producers to be added. This must be a list of producers. If multiple scopes are given,
        the producers are added to all scopes.
    
    """

    def add_producers(self, scopes, producers):
        if not isinstance(scopes, list):
            scopes = [scopes]
        if not isinstance(producers, list):
            producers = [producers]
        for scope in scopes:
            self.producers[scope].update(producers)
            self.available_outputs[scope].update(
                CollectProducersOutput(producers, scope)
            )

    """
    Function used to add outputs to the configuration. 

    Args:
        scopes: The scopes to which the outputs should be added. This can be a list of scopes or a single scope.
        outputs: The outputs to be added. This must be a list of outputs. If multiple scopes are given, 
        the outputs are added to all scopes.
    Returns:
        None
    
    """

    def add_outputs(self, scopes, output):
        if not isinstance(scopes, list):
            scopes = [scopes]
        if not isinstance(output, list):
            output = [output]
        for scope in scopes:
            self.outputs[scope].update(output)

    def add_shift(self, shift):
        if not isinstance(shift, SystematicShift):
            raise TypeError("shift must be of type SystematicShift")
        scopes_to_shift = [
            scope for scope in shift.get_scopes() if scope in self.scopes
        ]
        for scope in scopes_to_shift:
            shift.apply(scope)
            config_change = shift.get_shift_config(scope)
            self.shifts[scope][shift.name] = config_change

    def add_modification_rule(self, scopes, rule):
        if not isinstance(scopes, list):
            scopes = [scopes]
        rule.set_scopes(scopes)
        rule.set_global_scope(self.global_scope)
        self.rules.add(rule)

    def resolve_modifiers(self, parameters) -> dict:
        resolved = {}
        for key in parameters:
            resolved_value = None
            log.debug("Resolved: {}".format(resolved))
            log.debug("Parameter: {}".format(key))
            log.debug("Paramters: {}".format(parameters))
            if isinstance(parameters[key], dict):
                resolved.update(self.resolve_modifiers(parameters[key]))
            elif isinstance(parameters[key], SampleModifier):
                log.debug("Applying Samplemodifier")
                resolved_value = parameters[key].apply(self.sample)
            elif isinstance(parameters[key], EraModifier):
                log.debug("Applying Eramodifier")
                resolved_value = parameters[key].apply(self.era)
            else:
                resolved_value = parameters[key]
            resolved.update({key: resolved_value})
        return resolved

    def remove_empty_scopes(self):
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

    def optimize(self):
        self.apply_rules()
        self.remove_empty_scopes()
        for scope in self.scopes:
            log.debug("Optimizing Producer Ordering in scope {}".format(scope))
            ordering = ProducerOrdering(
                global_producers=self.producers[self.global_scope],
                scope=scope,
                producer_ordering=self.producers[scope],
            )
            ordering.Optimize()
            self.producers[scope] = ordering.ordering

    def validate_outputs(self):
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

    def apply_rules(self):
        for rule in self.rules:
            rule.apply(self.sample, self.producers, self.outputs)

    def validate(self):
        self.validate_outputs()

    def report(self):
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
        returnstr = "Configuration:"
        returnstr += "  Era: {}".format(self.era)
        returnstr += "  Sample: {}".format(self.sample)
        returnstr += "  Channels: {}".format(self.channels)
        returnstr += "  Shifts: {}".format(self.shifts)
        returnstr += "  Outputs:"
        for scope in self.scopes:
            returnstr += "    {}: {}".format(scope, self.outputs[scope])
        returnstr += "  Producers:"
        for scope in self.scopes:
            returnstr += "    {}: {}".format(scope, self.producers[scope])
        returnstr += "  Rules:"
        for scope in self.scopes:
            returnstr += "    {}: {}".format(scope, self.rules[scope])
        returnstr += "  Config Parameters:"
        for scope in self.scopes:
            returnstr += "    {}: {}".format(scope, self.config_parameters[scope])
        return ""

    def dump_dict(self) -> dict:
        returndict = {}
        returndict[""] = {}
        returndict["output"] = {}
        returndict["producers"] = {}
        for scope in self.scopes + [self.global_scope]:
            returndict[""][scope] = self.config_parameters[scope]
            returndict["producers"][scope] = self.producers[scope]
            if scope is not self.global_scope:
                log.info(sorted(list(self.outputs[scope])))
                returndict["output"][scope] = sorted(list(self.outputs[scope]))
                # returndict["output"][scope] = self.outputs[scope]
            # add systematic shifts
            for shift in self.shifts[scope]:
                log.warning("Adding shift {} in scope {}".format(shift, scope))
                log.warning("  {}".format(self.shifts[scope][shift]))
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


class Modifier(object):
    def __init__(self, modifier_dict):
        self.modifier_dict = modifier_dict

    def apply(self, configset):
        pass

    def __str__(self) -> str:
        return "Modifier: {}".format(self.modifier_dict)

    def __repr__(self) -> str:
        return "Modifier: {}".format(self.modifier_dict)


class SampleModifier(Modifier):
    def __init__(self, modifier_dict, default=None):
        super(SampleModifier, self).__init__(modifier_dict)
        self.modifier_dict = modifier_dict
        self.default = default
        self.samples = self.modifier_dict.keys()

    def apply(self, configsample):
        if configsample in self.samples:
            return self.modifier_dict[configsample]
        elif self.default is not None:
            return self.default
        else:
            raise SampleConfigurationError(configsample, self.samples)


class EraModifier(Modifier):
    def __init__(self, modifier_dict, default=None):
        super(EraModifier, self).__init__(modifier_dict)
        self.modifier_dict = modifier_dict
        self.default = default
        self.eras = self.modifier_dict.keys()

    def apply(self, configera):
        if configera in self.eras:
            return self.modifier_dict[configera]
        elif self.default is not None:
            return self.default
        else:
            raise EraConfigurationError(configera, self.eras)
