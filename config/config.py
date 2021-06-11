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
        "global": {
            "RunLumiEventFilter_Quantities": ["event"],
            "RunLumiEventFilter_Quantity_Types": ["ULong64_t"],
            "RunLumiEventFilter_Selections": ["271361"],
            "min_tau_pt": 30.0,
            "max_tau_eta": 2.3,
            "max_tau_dz": 0.2,
            "tau_dms": "0,1,10,11",
            "min_muon_pt": 10.0,
            "max_muon_eta": 2.4,
            "max_muon_dxy": 0.045,
            "max_muon_dz": 0.2,
            "min_jet_pt": 30,
            "max_jet_eta": 4.7,
            "jet_id": 2,  # second bit is tight JetID
            "min_bjet_pt": 20,
            "max_bjet_eta": 2.4,
            "btag_cut": 0.2783,  # medium
            "min_ele_pt": 10.0,
            "max_ele_eta": 2.4,
            "max_ele_dxy": 0.045,
            "max_ele_dz": 0.2,
            "max_ele_iso": 0.30,
            "ele_id": "Electron_mvaFall17V2noIso_WP90",
            "met_filters": ["Flag_goodVertices", "Flag_METFilters"],
            "tau_id": [
                "Tau_idDeepTau2017v2p1VSjet",
                "Tau_idDeepTau2017v2p1VSe",
                "Tau_idDeepTau2017v2p1VSmu",
            ],
            "tau_id_idx": [4, 4, 1],
            "muon_id": "Muon_mediumId",
            # "muon_iso": "Muon_pfRelIso04_all",
            "muon_iso_cut": 0.2,
            "JEC_shift_sources": '{""}',
            "JE_scale_shift": 0,
            "JE_reso_shift": 0,
            "tau_ES_shift_DM0": 1.0,
            "tau_ES_shift_DM1": 1.0,
            "tau_ES_shift_DM10": 1.0,
            "tau_ES_shift_DM11": 1.0,
        },
        "mt": {
            "mu_idx": 0,
            "min_muon_pt": 10.0,
            "max_muon_eta": 2.4,
            "muon_iso_cut": 0.15,
            "require_candidate": ["nTau", "nMuon"],
            "require_candidate_number": [1, 1],
            "deltaR_jet_veto": 0.5,
            "muon_sf_workspace": "data/muon_corrections/htt_scalefactors_legacy_2018_muons.root",
            "muon_sf_id_name": "m_id_kit_ratio",
            "muon_sf_id_args": "m_pt,m_eta",
            "muon_sf_iso_name": "m_iso_binned_kit_ratio",
            "muon_sf_iso_args": "m_pt,m_eta,m_iso",
        },
    }

    config = {"": base_config}

    config["producers"] = {
        "global": [
            # RunLumiEventFilter,
            Lumi,
            MetFilter,
            TauEnergyCorrection,
            GoodTaus,
            BaseMuons,
            BaseElectrons,
            JetEnergyCorrection,
            GoodJets,
            GoodBJets,
        ],
        "mt": [
            GoodMuons,
            MTPairSelection,
            GoodMTPairFilter,
            VetoMuons,
            ExtraMuonsVeto,
            ExtraElectronsVeto,
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
            q.muon_veto_flag,
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
        {"global": {"tau_ES_shift_DM0": 1.002}},
        [[TauPtCorrection, "global"]],
        [[LVMu1, "mt"], [VetoMuons, "mt"]],
    )
    AddSystematicShift(
        config,
        "tauES_1prong0pizeroDown",
        {"global": {"tau_ES_shift_DM0": 0.998}},
        [[TauPtCorrection, "global"]],
        [[LVMu1, "mt"], [VetoMuons, "mt"]],
    )
    # # Jet energy resolution
    # shift_dict = {"JE_reso_shift": 1}
    # AddSystematicShift(config, "jerUncUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_reso_shift": -1}
    # AddSystematicShift(config, "jerUncDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # # Jet energy scale
    # JEC_sources = '{"SinglePionECAL", "SinglePionHCAL", "AbsoluteMPFBias", "AbsoluteScale", "Fragmentation", "PileUpDataMC", "RelativeFSR", "PileUpPtRef"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncAbsoluteUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncAbsoluteDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])

    # JEC_sources = '{"AbsoluteStat", "TimePtEta", "RelativeStatFSR"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncAbsoluteYearUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]]
    # )
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncAbsoluteYearDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]]
    # )

    # JEC_sources = '{"FlavorQCD"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncFlavorQCDUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncFlavorQCDDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])

    # JEC_sources = '{"PileUpPtEC1", "PileUpPtBB", "RelativePtBB"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncBBEC1Up", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncBBEC1Down", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])

    # JEC_sources = '{"RelativeJEREC1", "RelativePtEC1", "RelativeStatEC"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncBBEC1YearUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncBBEC1YearDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])

    # JEC_sources = '{"RelativePtHF", "PileUpPtHF", "RelativeJERHF"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncHFUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncHFDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])

    # JEC_sources = '{"RelativeStatHF"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncHFYearUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncHFYearDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])

    # JEC_sources = '{"PileUpPtEC2"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncEC2Up", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncEC2Down", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])

    # JEC_sources = '{"RelativeJEREC2", "RelativePtEC2"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUnjecUncEC2YearUpcHFUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]]
    # )
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncEC2YearDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])

    # JEC_sources = '{"RelativeBal"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(config, "jecUncRelativeBalUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]])
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncRelativeBalDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]]
    # )

    # JEC_sources = '{"RelativeSample"}'
    # shift_dict = {"JE_scale_shift": 1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncRelativeSampleYearUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]]
    # )
    # shift_dict = {"JE_scale_shift": -1, "JEC_shift_sources": JEC_sources}
    # AddSystematicShift(
    #     config, "jecUncRelativeSampleYearDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]]
    # )

    return config
