import copy


def AddSystematicShift(
    config, name, change_dict, base_producers, sanetize_producers=[], scope="global"
):
    name = "__" + name
    shift_config = copy.deepcopy(config[""])
    shift_config.update(change_dict)
    config[name] = shift_config
    for producer in sanetize_producers:
        producer[0].ignore_shift(name, producer[1])
    for producer in base_producers:
        producer.shift(name, scope)
