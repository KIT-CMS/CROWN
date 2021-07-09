from code_generation.producers.jets import *
from code_generation.producers.taus import *
from code_generation.producers.muons import *
from code_generation.producers.electrons import *
from code_generation.producers.genparticles import *
from code_generation.producers.pairselection import *
from code_generation.producers.pairquantities import *
from code_generation.producers.event import *
from code_generation.producers.scalefactors import *
from code_generation.producers.triggers import *
import code_generation.quantities.output as q
from config.utility import (
    AddSystematicShift,
    ResolveSampleDependencies,
    ResolveEraDependencies,
    RemoveProducer,
    AppendProducer,
)


def build_config(era, sample):
    base_config = {
        "global": {
            "RunLumiEventFilter_Quantities": ["event"],
            "RunLumiEventFilter_Quantity_Types": ["ULong64_t"],
            "RunLumiEventFilter_Selections": ["271361"],
            "PU_reweighting_file": {
                "ERA_2016": "data/pileup/Data_Pileup_2016_271036-284044_13TeVMoriond17_23Sep2016ReReco_69p2mbMinBiasXS.root",
                "ERA_2017": "data/pileup/Data_Pileup_2017_294927-306462_13TeVSummer17_PromptReco_69p2mbMinBiasXS.root",
                "ERA_2018": "data/pileup/Data_Pileup_2018_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18.root",
            },
            "PU_reweighting_hist": "pileup",
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
            "max_ele_eta": 2.5,
            "max_ele_dxy": 0.045,
            "max_ele_dz": 0.2,
            "max_ele_iso": 0.3,
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
            "muon_iso_cut": 0.3,
            "JEC_shift_sources": '{""}',
            "JE_scale_shift": 0,
            "JE_reso_shift": 0,
            "tau_ES_shift_DM0": 1.0,
            "tau_ES_shift_DM1": 1.0,
            "tau_ES_shift_DM10": 1.0,
            "tau_ES_shift_DM11": 1.0,
            "min_dielectronveto_pt": 15.0,
            "dielectronveto_id": "Electron_cutBased",
            "dielectronveto_id_wp": 1,
            "min_dimuonveto_pt": 15.0,
            "dimuonveto_id": "Muon_looseId",
            "dileptonveto_dR": 0.15,
        },
        "mt": {
            "singlemoun_trigger": {
                "ERA_2018": [
                    {
                        "flagname": "singlemuon_24",
                        "hlt_path": "HLT_IsoMu24",
                        "ptcut": 25,
                        "etacut": 2.5,
                        "filterbit": 4,
                        "trigger_particle_id": 13,
                        "max_deltaR_triggermatch": 0.4,
                    },
                    {
                        "flagname": "singlemuon_27",
                        "hlt_path": "HLT_IsoMu27",
                        "ptcut": 28,
                        "etacut": 2.5,
                        "filterbit": 4,
                        "trigger_particle_id": 13,
                        "max_deltaR_triggermatch": 0.4,
                    },
                ],
            },
            "cross_trigger": {
                "ERA_2018": [
                    {
                        "flagname": "trg_crossmuon_mu20tau27_hps",
                        "hlt_path": "HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1",
                        "p1_ptcut": 21,
                        "p2_ptcut": 32,
                        "p1_etacut": 2.5,
                        "p2_etacut": 2.1,
                        "p1_filterbit": 4,
                        "p1_trigger_particle_id": 13,
                        "p2_filterbit": 0,
                        "p2_trigger_particle_id": 15,
                        "max_deltaR_triggermatch": 0.4,
                    }
                ],
            },
            "mu_idx": 0,
            "min_muon_pt": 23.0,
            "max_muon_eta": 2.1,
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
    all_channels = {
        "ggHNNLOweightsRootfile": "data/htxs/NNLOPS_reweight.root",
        "ggH_generator": "powheg",
    }
    for channel in ["mt"]:  # add em et tt here as soon as they appear in config
        base_config[channel].update(all_channels)

    config = {"": base_config}

    config["producers"] = {
        "global": [
            # RunLumiEventFilter,
            Lumi,
            MetFilter,
            PUweights,
            TauEnergyCorrection,
            GoodTaus,
            BaseMuons,
            BaseElectrons,
            DiLeptonVeto,
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
            GenerateSingleMuonTriggerFlags,
            GenerateCrossTriggerFlags,
        ],
    }

    config["producer_modifiers"] = [
        RemoveProducer(producers=[MuonIDIso_SF], samples=["data"], scopes=["mt"]),
        RemoveProducer(
            producers=[PUweights], samples=["data", "emb"], scopes=["global"]
        ),
        AppendProducer(
            producers=[GGH_NNLO_Reweighting], samples=["ggh"], scopes=["mt"]
        ),
        AppendProducer(
            producers=[GGH_WG1_Uncertainties], samples=["ggh"], scopes=["mt"]
        ),
        AppendProducer(
            producers=[QQH_WG1_Uncertainties], samples=["qqh"], scopes=["mt"]
        ),
        AppendProducer(producers=[TopPtReweighting], samples=["ttbar"], scopes=["mt"]),
    ]

    config["output"] = {
        "mt": [
            nanoAOD.run,
            q.lumi,
            nanoAOD.event,
            q.puweight,
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
            q.dimuon_veto,
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
            GenerateSingleMuonTriggerFlags.output_group,
            GenerateCrossTriggerFlags.output_group,
            nanoAOD.HTXS_Higgs_pt,
            nanoAOD.HTXS_njets30,
            nanoAOD.HTXS_stage_0,
            nanoAOD.HTXS_stage1_2_cat_pTjet30GeV,
        ]
    }

    for modifier in config["producer_modifiers"]:
        modifier.apply(sample, config["producers"], config["output"])
    ResolveSampleDependencies(config, sample)
    ResolveEraDependencies(config, era)

    AddSystematicShift(
        config,
        "tauES_1prong0pizeroUp",
        {"global": {"tau_ES_shift_DM0": 1.002}},
        [[TauPtCorrection, "global"]],
        sanetize_producers=[[LVMu1, "mt"], [VetoMuons, "mt"]],
    )
    AddSystematicShift(
        config,
        "tauES_1prong0pizeroDown",
        {"global": {"tau_ES_shift_DM0": 0.998}},
        [[TauPtCorrection, "global"]],
        sanetize_producers=[[LVMu1, "mt"], [VetoMuons, "mt"]],
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
