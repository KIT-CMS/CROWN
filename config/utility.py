import copy
import logging

log = logging.getLogger(__name__)


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


class ProducerRule:
    def __init__(self, producers, samples, scopes, invert=False, update_output=True):
        self.producers = producers
        self.samples = samples
        self.scopes = scopes
        self.invert = invert
        self.update_output = update_output

    def applicable(self, sample):
        applicable = sample in self.samples
        if self.invert:
            applicable = not applicable
        return applicable

    def operate(self, item, item_list):
        log.error("Operation not implemented for ProducerRule base class!")
        raise Exception

    def apply(self, sample, producer_dict, output_dict=None):
        if self.applicable(sample):
            for scope in self.scopes:
                for prod in self.producers:
                    self.operate(prod, producer_dict[scope])
                    if self.update_output and output_dict != None:
                        for q in CollectProducerOutput(prod):
                            self.operate(q, output_dict[scope])


class Prod_Remove(ProducerRule):
    def operate(self, item, item_list):
        item_list.remove(item)


class Prod_Append(ProducerRule):
    def operate(self, item, item_list):
        item_list.append(item)
