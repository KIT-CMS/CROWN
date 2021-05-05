import copy


def AddSystematicShift(config, name, change_dict, base_producers):
    name = "__" + name
    shift_config = copy.deepcopy(config[""])
    shift_config.update(change_dict)
    config[name] = shift_config
    for producer in base_producers:
        producer.shift(name)
