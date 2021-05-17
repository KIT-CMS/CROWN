import code_generation.producer as p
from code_generation.producers.jets import *
from code_generation.producer import q
from config.utility import AddSystematicShift


def build_config():
    base_config = {
        "min_tau_pt": 30.0,
        "max_tau_eta": 2.3,
        "max_tau_dz": 0.2,
        "min_muon_pt": 23.0,
        "max_muon_eta": 2.5,
        "min_jet_pt": 30,
        "max_jet_eta": 4.7,
        "jet_id": 2,  # second bit is tight JetID
        "deltaR_jet_veto": 0.5,
        "min_VetoElectron_pt": 10.0,
        "max_VetoElectron_eta": 2.5,
        "max_VetoElectron_iso": 0.30,
        "VetoElectron_id": "Electron_mvaFall17V2noIso_WP90",
        "met_filters": ["Flag_goodVertices", "Flag_METFilters"],
        "tau_id": [
            "Tau_idDeepTau2017v2p1VSjet",
            "Tau_idDeepTau2017v2p1VSe",
            "Tau_idDeepTau2017v2p1VSmu",
        ],
        "tau_id_idx": [4, 4, 1],
        "muon_id": "Muon_mediumId",
        # "muon_iso": "Muon_pfRelIso04_all",
        "muon_iso_cut": 0.15,
        "require_candidate": ["nTau", "nMuon"],
        "require_candidate_number": [1, 1],
    }

    config = {"": base_config}

    config["producers"] = {
        "global": [
            p.MetFilter,
            p.GoodTaus,
            p.GoodMuons,
            p.GoodElectronsVeto,
            # GoodJets,
        ],
        "mt": [
            p.MTPairSelection,
            p.GoodMTPairFilter,
            p.LVMu1,
            p.LVTau2,
            p.DiTauPairQuantities,
            VetoJets,
            BasicJetQuantities,
        ],
    }

    config["output"] = {
        "mt": [
            q.pt_1,
            q.pt_2,
            q.njets,
            q.jpt_1,
            q.jpt_2,
            q.jeta_1,
            q.jeta_2,
            q.jphi_1,
            q.jphi_2,
            q.mjj,
            q.electron_veto_flag,
            q.good_jet_collection,
        ]
    }

    shift_dict = {"min_tau_pt": 31.0}
    AddSystematicShift(config, "tauCutUp", shift_dict, [p.TauPtCut])

    return config
