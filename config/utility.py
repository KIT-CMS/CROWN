import copy


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
