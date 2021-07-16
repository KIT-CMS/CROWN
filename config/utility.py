from code_generation.optimizer import ProducerOrdering
import copy
import logging

log = logging.getLogger(__name__)

# Function for introducing systematic variations to producers and depending quantities
def AddSystematicShift(
    config, name, change_dict, base_producers, sanetize_producers=[]
):
    name = "__" + name
    shift_config = copy.deepcopy(config[""])
    for scope in change_dict:
        shift_config[scope].update(change_dict[scope])
    config[name] = shift_config
    for producer in sanetize_producers:
        producer[0].ignore_shift(name, producer[1])
    for producer in base_producers:
        producer[0].shift(name, producer[1])


# Function for introducing systematic variations to producers and depending quantities by adding an already shifted input quantity
def SystematicShiftByInputQuantity(
    config, shiftname, base_producers, external_dict
):
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
                else:
                    log.error("Era {} not available in {}!".format(era, key))
                    raise Exception
            else:
                ResolveEraDependencies(config[key], era)


# Helper function to retrieve all quantities produced by a producer or producer group
def CollectProducerOutput(producer):
    output = []
    if producer.output != None:
        if isinstance(producer.output, list):
            output.extend(producer.output)
        else:
            output.append(producer.output)
    if hasattr(producer, "producers"):
        for prod in producer.producers:
            output.extend(CollectProducerOutput(prod))
    return output


# Base class of modifiers that are used to apply sample specific modifications to the producer lists
class ProducerRule:
    def __init__(self, producers, samples, scopes, invert=False, update_output=True):
        self.producers = producers
        self.samples = samples
        self.scopes = scopes
        self.invert = invert
        self.update_output = update_output

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
        if self.applicable(sample):
            for scope in self.scopes:
                for prod in self.producers:
                    self.operate(prod, producer_dict[scope])
                    if self.update_output and output_dict != None:
                        for q in CollectProducerOutput(prod):
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
        item_list.remove(item)


# Modifier class that can append producers to lists
class AppendProducer(ProducerRule):
    def operate(self, item, item_list):
        item_list.append(item)


def OptimizeProducerOrdering(config):
    scopes = config["producers"].keys()
    for scope in scopes:
        log.info("Optimizing Producer Ordering in scope {}".format(scope))
        ordering = ProducerOrdering(config, scope)
        ordering.Optimize()
        config["producers"][scope] = ordering.ordering
