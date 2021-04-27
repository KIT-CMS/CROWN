import copy

import code_generation.producer as p


def build_config():
    base_config = {
        "ptcut": 2.0,
        "etacut": 1.0,
        "met_filters": ["Flag_goodVertices", "Flag_METFilters"],
    }

    config = {"": copy.deepcopy(base_config), "_tauEsUp": copy.deepcopy(base_config)}

    config["producers"] = ["MetFilter"]

    config["_tauEsUp"]["ptcut"] = 2.1
    # write some modifier tools for creating shifts. Should these automatically determine producers that consume params and shift these?
    # Or simply add an explicit command to the modifier tool to shift a producer
    # AddShift(config, "_tauEsDown", dict)

    return config
