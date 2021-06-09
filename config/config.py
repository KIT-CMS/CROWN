from code_generation.filters.filters import *
from code_generation.producers.jets import *
from code_generation.producers.taus import *
from code_generation.producers.muons import *
from code_generation.producers.electrons import *
from code_generation.producers.genparticles import *
from code_generation.producers.pairselection import *
from code_generation.producers.pairquantities import *
from code_generation.producers.event import *
from code_generation.producers.scalefactors import *
import code_generation.quantities.output as q
from config.utility import AddSystematicShift


def build_config():
    base_config = {
        "min_tau_pt": 30.0,
        "max_tau_eta": 2.3,
        "max_tau_dz": 0.2,
        "tau_dms": "0,1,10,11",
        "min_muon_pt": 23.0,
        "max_muon_eta": 2.5,
        "max_muon_dxy": 0.045,
        "max_muon_dz": 0.2,
        "min_jet_pt": 30,
        "max_jet_eta": 4.7,
        "jet_id": 2,  # second bit is tight JetID
        "min_bjet_pt": 20,
        "max_bjet_eta": 2.4,
        "btag_cut": 0.2783,  # medium
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
        "JEC_shift_sources": '{""}',
        "JE_scale_shift": 0,
        "JE_reso_shift": 0,
        "tau_ES_shift_DM0": 1.0,
        "tau_ES_shift_DM1": 1.0,
        "tau_ES_shift_DM10": 1.0,
        "tau_ES_shift_DM11": 1.0,
        "muon_sf_workspace": "data/muon_corrections/htt_scalefactors_legacy_2018_muons.root",
        "muon_sf_id_name": "m_id_kit_ratio",
        "muon_sf_id_args": "m_pt,m_eta",
        "muon_sf_iso_name": "m_iso_binned_kit_ratio",
        "muon_sf_iso_args": "m_pt,m_eta,m_iso",
        "RunLumiEventFilter_Quantities": ["event"],
        "RunLumiEventFilter_Quantity_Types": ["ULong64_t"],
        "RunLumiEventFilter_Selections": ["271361"],
    }

    config = {"": base_config}

    config["producers"] = {
        "global": [
            # RunLumiEventFilter,
            Lumi,
            MetFilter,
            TauEnergyCorrection,
            GoodTaus,
            GoodMuons,
            GoodElectronsVeto,
            JetEnergyCorrection,
            GoodJets,
            GoodBJets,
        ],
        "mt": [
            MTPairSelection,
            GoodMTPairFilter,
            LVMu1,
            LVTau2,
            DiTauPairQuantities,
            JetCollection,
            BasicJetQuantities,
            BJetCollection,
            BasicBJetQuantities,
            GenDiTauPairQuantities,
            MuonIDIso_SF,
        ],
    }

    config["output"] = {
        "mt": [
            nanoAOD.run,
            q.lumi,
            nanoAOD.event,
            q.pt_1,
            q.pt_2,
            q.eta_1,
            q.eta_2,
            q.phi_1,
            q.phi_2,
            q.njets,
            q.jpt_1,
            q.jpt_2,
            q.jeta_1,
            q.jeta_2,
            q.jphi_1,
            q.jphi_2,
            q.mjj,
            q.m_vis,
            q.electron_veto_flag,
            q.nbtag,
            q.bpt_1,
            q.bpt_2,
            q.beta_1,
            q.beta_2,
            q.bphi_1,
            q.bphi_2,
            q.mass_1,
            q.mass_2,
            q.dxy_1,
            q.dxy_2,
            q.dz_1,
            q.dz_2,
            q.q_1,
            q.q_2,
            q.iso_1,
            q.iso_2,
            q.decaymode_2,
            q.gen_match_2,
            q.gen_pt_1,
            q.gen_eta_1,
            q.gen_phi_1,
            q.gen_mass_1,
            q.gen_pdgid_1,
            q.gen_pt_2,
            q.gen_eta_2,
            q.gen_phi_2,
            q.gen_mass_2,
            q.gen_pdgid_2,
            q.gen_m_vis,
            q.taujet_pt_2,
            q.gen_taujet_pt_2,
            q.idWeight_1,
            q.isoWeight_1,
        ]
    }

    AddSystematicShift(
        config,
        "tauES_1prong0pizeroUp",
        {"tau_ES_shift_DM0": 1.002},
        [TauPtCorrection],
        [[LVMu1, "mt"]],
    )
    AddSystematicShift(
        config,
        "tauES_1prong0pizeroDown",
        {"tau_ES_shift_DM0": 0.998},
        [TauPtCorrection],
        [[LVMu1, "mt"]],
    )
    # # Jet energy resolution
    # shift_dict = {"JE_reso_shift": 1}
    # AddSystematicShift(config, "jerUncUp", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_reso_shift": -1}
    # AddSystematicShift(config, "jerUncDown", shift_dict, [JetEnergyCorrection])
    # # Jet energy scale
    # JEC_sources = '{"SinglePionECAL", "SinglePionHCAL", "AbsoluteMPFBias", "AbsoluteScale", "Fragmentation", "PileUpDataMC", "RelativeFSR", "PileUpPtRef"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncAbsoluteUp", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncAbsoluteDown", shift_dict, [JetEnergyCorrection])

    # JEC_sources = '{"AbsoluteStat", "TimePtEta", "RelativeStatFSR"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncAbsoluteYearUp", shift_dict, [JetEnergyCorrection]
    # )
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncAbsoluteYearDown", shift_dict, [JetEnergyCorrection]
    # )

    # JEC_sources = '{"FlavorQCD"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncFlavorQCDUp", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncFlavorQCDDown", shift_dict, [JetEnergyCorrection])

    # JEC_sources = '{"PileUpPtEC1", "PileUpPtBB", "RelativePtBB"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncBBEC1Up", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncBBEC1Down", shift_dict, [JetEnergyCorrection])

    # JEC_sources = '{"RelativeJEREC1", "RelativePtEC1", "RelativeStatEC"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncBBEC1YearUp", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncBBEC1YearDown", shift_dict, [JetEnergyCorrection])

    # JEC_sources = '{"RelativePtHF", "PileUpPtHF", "RelativeJERHF"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncHFUp", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncHFDown", shift_dict, [JetEnergyCorrection])

    # JEC_sources = '{"RelativeStatHF"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncHFYearUp", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncHFYearDown", shift_dict, [JetEnergyCorrection])

    # JEC_sources = '{"PileUpPtEC2"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncEC2Up", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncEC2Down", shift_dict, [JetEnergyCorrection])

    # JEC_sources = '{"RelativeJEREC2", "RelativePtEC2"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUnjecUncEC2YearUpcHFUp", shift_dict, [JetEnergyCorrection]
    # )
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncEC2YearDown", shift_dict, [JetEnergyCorrection])

    # JEC_sources = '{"RelativeBal"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncRelativeBalUp", shift_dict, [JetEnergyCorrection])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncRelativeBalDown", shift_dict, [JetEnergyCorrection]
    # )

    # JEC_sources = '{"RelativeSample"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncRelativeSampleYearUp", shift_dict, [JetEnergyCorrection]
    # )
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncRelativeSampleYearDown", shift_dict, [JetEnergyCorrection]
    # )

    return config
