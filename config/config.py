import copy

import code_generation.producer as p


def build_config():
    base_config = {
        "ptcut": 30.0,
        "etacut": 2.3,
        "dzcut": 0.2,
        "met_filters": ["Flag_goodVertices", "Flag_METFilters"],
        "tau_id": ["Tau_idDeepTau2017v2p1VSjet", "Tau_idDeepTau2017v2p1VSe", "Tau_idDeepTau2017v2p1VSmu"],
        "tau_id_idx": [4, 4, 1]
    }

    config = {"": copy.deepcopy(base_config), "_tauCutUp": copy.deepcopy(base_config)}

    config["producers"] = ["MetFilter", "GoodTaus"]

    config["_tauCutUp"]["ptcut"] = 2.1
    config["_tauCutUp"]["shiftbase"] = ["TauPtCut"]
    # write some modifier tools for creating shifts. Should these automatically determine producers that consume params and shift these?
    # Or simply add an explicit command to the modifier tool to shift a producer
    # AddShift(config, "_tauEsDown", dict)

    return config
