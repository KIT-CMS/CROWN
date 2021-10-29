from code_generation.producer import ProducerGroup
from code_generation.optimizer import ProducerOrdering
from code_generation.quantity import NanoAODQuantity, Quantity
from code_generation.exceptions import *
import copy
import logging


log = logging.getLogger(__name__)

# Function for introducing systematic variations to producers and depending quantities
def AddSystematicShift(
    config, name, change_dict, base_producers, sanitize_producers=[]
):
    name = "__" + name
    shift_config = copy.deepcopy(config[""])
    for scope in change_dict:
        shift_config[scope].update(change_dict[scope])
    config[name] = shift_config
    for producer in sanitize_producers:
        producer[0].ignore_shift(name, producer[1])
    for producer in base_producers:
        producer[0].shift(name, producer[1])


# Function for introducing systematic variations to producers and depending quantities by adding an already shifted input quantity
def SystematicShiftByInputQuantity(config, shiftname, external_dict):
    shiftname = "__" + shiftname
    config[shiftname] = copy.deepcopy(config[""])
    for quantity in external_dict.keys():
        quantity.register_external_shift(
            shift_name=shiftname,
            external_name=external_dict[quantity],
        )


# Function for resolving sample dependencies in the config dict for a given sample
def ResolveSampleDependencies(config, sample):
    for key in config.keys():
        if isinstance(config[key], dict):
            if list(config[key].keys())[0].startswith("SAMPLE_"):
                if "SAMPLE_" + sample in config[key].keys():
                    config[key] = config[key]["SAMPLE_" + sample]
                elif "SAMPLE_DEFAULT" in config[key].keys():
                    config[key] = config[key]["SAMPLE_DEFAULT"]
                else:
                    log.error("Sample {} not available in {}!".format(sample, key))
                    raise Exception
            else:
                ResolveSampleDependencies(config[key], sample)


# Function for resolving era dependencies in the config dict for a given era
def ResolveEraDependencies(config, era):
    for key in config.keys():
        if isinstance(config[key], dict):
            if list(config[key].keys())[0].startswith("ERA_"):
                if "ERA_" + era in config[key].keys():
                    config[key] = config[key]["ERA_" + era]
                elif "ERA_DEFAULT" in config[key].keys():
                    config[key] = config[key]["ERA_DEFAULT"]
                else:
                    log.error("Era {} not available in {}!".format(era, key))
                    raise Exception
            else:
                ResolveEraDependencies(config[key], era)


# Helper function to retrieve all quantities produced by a producer or producer group
def CollectProducerOutput(producer, scope):
    output = []
    if producer.output is not None:
        if isinstance(producer.output, list):
            output.extend(producer.output)
        else:
            output.append(producer.output)
    if hasattr(producer, "producers"):
        for prod in producer.producers[scope]:
            output.extend(CollectProducerOutput(prod, scope))
    return output


def CollectProducersOutput(producers, scope):
    output = []
    for producer in producers:
        if producer.output is not None:
            if isinstance(producer.output, list):
                output.extend(producer.output)
            else:
                output.append(producer.output)
        if hasattr(producer, "producers"):
            for prod in producer.producers[scope]:
                output.extend(CollectProducerOutput(prod, scope))
    return output


def ExpandProducerConfig(producers, scope):
    # we expand all producers groups to individual producers
    full_producerlist = []
    for producer in producers:
        if isinstance(producer, ProducerGroup):
            full_producerlist.extend(
                ExpandProducerConfig(producer.producers[scope], scope)
            )
        else:
            full_producerlist.append(producer)
    # log.info(full_producerlist)
    return full_producerlist


# Base class of modifiers that are used to apply sample specific modifications to the producer lists
class ProducerRule:
    def __init__(self, producers, samples, scopes, invert=False, update_output=True):
        self.producers = producers
        self.samples = samples
        self.scopes = scopes
        self.invert = invert
        self.update_output = update_output

    def __str__(self) -> str:
        return "ProducerRule - update {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    def __repr__(self) -> str:
        return "ProducerRule - update {} for {} in scopes {}".format(
            self.producers, self.samples, self.scopes
        )

    # Evaluate whether modification should be applied depending on sample and inversion flag
    def applicable(self, sample):
        applicable = sample in self.samples
        if self.invert:
            applicable = not applicable
        return applicable

    # Placeholder for the actual operation on a list. To be overwritten by inheriting classes
    def operate(self, item, item_list):
        log.error("Operation not implemented for ProducerRule base class!")
        raise Exception

    # Method to be called at the level of config generation in order to evaluate conditions and apply the modification
    def apply(self, sample, producer_dict, output_dict=None):
        log.debug("Applying {} for {}".format(self, sample))
        if self.applicable(sample):
            for scope in self.scopes:
                for prod in self.producers:
                    self.operate(prod, producer_dict[scope])
                for prod in self.producers:
                    if self.update_output and output_dict is not None:
                        for q in CollectProducerOutput(prod, scope):
                            scopelist = (
                                [scope]
                                if scope in output_dict.keys()
                                else output_dict.keys()
                                if scope == "global"
                                else []
                            )
                            for scope2 in scopelist:
                                self.operate(q, output_dict[scope2])


# Modifier class that can remove producers from lists
class RemoveProducer(ProducerRule):
    def operate(self, item, item_list):
        try:
            log.debug(
                "RemoveProducer: Removing {} from list {}".format(item, item_list)
            )
            item_list.remove(item)
        except ValueError:
            log.warning("Error when applying {} ".format(self))
            log.warning("Cannot remove {} from {}".format(item, item_list))
            pass


# Modifier class that can append producers to lists
class AppendProducer(ProducerRule):
    def operate(self, item, item_list):
        log.debug("AppendProducer: Adding {} to list {}".format(item, item_list))
        item_list.append(item)


def OptimizeProducerOrdering(config):
    scopes = config["producers"].keys()
    for scope in scopes:
        log.info("Optimizing Producer Ordering in scope {}".format(scope))
        ordering = ProducerOrdering(config, scope)
        ordering.Optimize()
        config["producers"][scope] = ordering.ordering


"""
This function is used to validate the set of outputs defined in the config.
If the output list contains quantities that are not produced by any producer,
the function will raise an exception.

Args:
    config: The config dict
"""


def ValidateOutputs(config):
    global_producers = ExpandProducerConfig(config["producers"]["global"], "global")
    for scope in set(config["producers"]) - {"global"}:
        producers = global_producers + ExpandProducerConfig(
            config["producers"][scope], scope
        )
        provided_outputs = [
            producer.output for producer in producers if producer.output is not None
        ]
        provided_outputs = set(
            output
            for sublist in provided_outputs
            for output in sublist
            if not isinstance(output, NanoAODQuantity)
        )
        required_outputs = set(
            output
            for output in config["output"][scope]
            if not isinstance(output, NanoAODQuantity)
        )
        # find all outputs that are required but not provided
        missing_outputs = required_outputs - provided_outputs
        if len(missing_outputs) > 0:
            raise InvalidOutputError(scope, missing_outputs)


def SetSampleParameters(config, sample, available_sample_types, channels):
    sample_parameters = {}
    for sampletype in available_sample_types:
        if sample == sampletype:
            sample_parameters["is_{}".format(sampletype)] = True
        else:
            sample_parameters["is_{}".format(sampletype)] = False
    for channel in channels:
        config[channel].update(sample_parameters)


def SetCommonParameters(base_config, commons, channels):
    for channel in channels:
        base_config[channel].update(commons)
    return {"": base_config}
