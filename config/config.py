import copy

import code_generation.producer as p


def build_config():
    base_config = {
        "min_tau_pt": 30.0,
        "max_tau_eta": 2.3,
        "max_tau_dz": 0.2,
        "min_muon_pt": 23.0,
        "max_muon_eta": 2.5,
        "met_filters": ["Flag_goodVertices", "Flag_METFilters"],
        "tau_id": [
            "Tau_idDeepTau2017v2p1VSjet",
            "Tau_idDeepTau2017v2p1VSe",
            "Tau_idDeepTau2017v2p1VSmu",
        ],
        "tau_id_idx": [4, 4, 1],
        "muon_id": "Muon_mediumId",
        "muon_iso": "Muon_pfRelIso04_all",
        "muon_iso_cut": 0.15,
        "require_candidate": ["nTau", "nMuon"],
        "require_candidate_number": [1, 1],
    }

    config = {"": copy.deepcopy(base_config), "_tauCutUp": copy.deepcopy(base_config)}

    config["producers"] = ["MetFilter", "GoodTaus", "GoodMuons"]

    config["_tauCutUp"]["min_tau_pt"] = 31.0
    config["_tauCutUp"]["shiftbase"] = ["TauPtCut"]
    # write some modifier tools for creating shifts. Should these automatically determine producers that consume params and shift these?
    # Or simply add an explicit command to the modifier tool to shift a producer
    # AddShift(config, "_tauEsDown", dict)

    return config
