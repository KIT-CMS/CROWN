from code_generation.configuration import (
    Configuration,
    SampleModifier,
    EraModifier,
)
from code_generation.rules import AppendProducer, RemoveProducer
from code_generation.systematics import SystematicShift
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
from code_generation.producers.met import *
import code_generation.quantities.output as q


def build_config(
    era,
    sample,
    channels,
    shifts,
    available_sample_types,
    available_eras,
    available_channels,
):

    configuration = Configuration(
        era,
        sample,
        channels,
        shifts,
        available_sample_types,
        available_eras,
        available_channels,
    )
    # first add default parameters necessary for all scopes
    configuration.add_config_parameters(
        "global",
        {
            "RunLumiEventFilter_Quantities": ["event"],
            "RunLumiEventFilter_Quantity_Types": ["ULong64_t"],
            "RunLumiEventFilter_Selections": ["271361"],
            "PU_reweighting_file": EraModifier(
                {
                    "2016": "data/pileup/Data_Pileup_2016_271036-284044_13TeVMoriond17_23Sep2016ReReco_69p2mbMinBiasXS.root",
                    "2017": "data/pileup/Data_Pileup_2017_294927-306462_13TeVSummer17_PromptReco_69p2mbMinBiasXS.root",
                    "2018": "data/pileup/Data_Pileup_2018_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18.root",
                }
            ),
            "golden_json_file": EraModifier(
                {
                    "2016": "data/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
                    "2017": "data/golden_json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
                    "2018": "data/golden_json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
                }
            ),
            "PU_reweighting_hist": "pileup",
            "met_filters": ["Flag_goodVertices", "Flag_METFilters"],
        },
    )
    # Tau base selection:
    configuration.add_config_parameters(
        "global",
        {
            "min_tau_pt": 30.0,
            "max_tau_eta": 2.3,
            "max_tau_dz": 0.2,
            "tau_dms": "0,1,10,11",
            # "tau_id": [
            #     "Tau_idDeepTau2017v2p1VSjet",
            #     "Tau_idDeepTau2017v2p1VSe",
            #     "Tau_idDeepTau2017v2p1VSmu",
            # ],
            # "tau_id_idx": [4, 4, 1],
            "vsjet_tau_id_bit": 4,
            "vsele_tau_id_bit": 4,
            "vsmu_tau_id_bit": 1,
            "tau_ES_shift_DM0": 1.0,
            "tau_ES_shift_DM1": 1.0,
            "tau_ES_shift_DM10": 1.0,
            "tau_ES_shift_DM11": 1.0,
        },
    )
    # muon base selection:
    configuration.add_config_parameters(
        "global",
        {
            "min_muon_pt": 10.0,
            "max_muon_eta": 2.4,
            "max_muon_dxy": 0.045,
            "max_muon_dz": 0.2,
            "muon_id": "Muon_mediumId",
            "muon_iso_cut": 0.3,
        },
    )
    # electron base selection:
    configuration.add_config_parameters(
        "global",
        {
            "min_ele_pt": 10.0,
            "max_ele_eta": 2.5,
            "max_ele_dxy": 0.045,
            "max_ele_dz": 0.2,
            "max_ele_iso": 0.3,
            "ele_id": "Electron_mvaFall17V2noIso_WP90",
        },
    )
    # jet base selection:
    configuration.add_config_parameters(
        "global",
        {
            "min_jet_pt": 30,
            "max_jet_eta": 4.7,
            "jet_id": 2,  # second bit is tight JetID
            "JEC_shift_sources": '{""}',
            "JE_scale_shift": 0,
            "JE_reso_shift": 0,
        },
    )
    # bjet base selection:
    configuration.add_config_parameters(
        "global",
        {
            "min_bjet_pt": 20,
            "max_bjet_eta": 2.4,
            "btag_cut": 0.2783,  # medium
        },
    )
    # leptonveto base selection:
    configuration.add_config_parameters(
        "global",
        {
            "min_dielectronveto_pt": 15.0,
            "dielectronveto_id": "Electron_cutBased",
            "dielectronveto_id_wp": 1,
            "min_dimuonveto_pt": 15.0,
            "dimuonveto_id": "Muon_looseId",
            "dileptonveto_dR": 0.15,
        },
    )
    ###### Channel Specifics ######
    ## MT/MM channel Muon selection
    configuration.add_config_parameters(
        ["mt", "mm"],
        {
            "mu_idx": 0,
            "min_muon_pt": 23.0,
            "max_muon_eta": 2.1,
            "muon_iso_cut": 0.15,
            "muon_sf_workspace": "data/muon_corrections/htt_scalefactors_legacy_2018_muons.root",
            "muon_sf_id_name": "m_id_kit_ratio",
            "muon_sf_id_args": "m_pt,m_eta",
            "muon_sf_iso_name": "m_iso_binned_kit_ratio",
            "muon_sf_iso_args": "m_pt,m_eta,m_iso",
        },
    )
    ## MT/MM channel misc settings
    configuration.add_config_parameters(
        ["mt", "mm"],
        {
            "deltaR_jet_veto": 0.5,
        },
    )
    ## MT/MM channel MET selection
    configuration.add_config_parameters(
        ["mt", "mm"],
        {
            "propagateLeptons": SampleModifier(
                {"data": False, "emb": False},
                default=True,
            ),
            "propagateJets": SampleModifier(
                {"data": False, "emb": False},
                default=True,
            ),
            "recoil_corrections_file": EraModifier(
                {
                    "2016": "data/recoil_corrections/Type1_PuppiMET_2016.root",
                    "2017": "data/recoil_corrections/Type1_PuppiMET_2017.root",
                    "2018": "data/recoil_corrections/Type1_PuppiMET_2018.root",
                }
            ),
            "recoil_systematics_file": EraModifier(
                {
                    "2016": "data/recoil_corrections/PuppiMETSys_2016.root",
                    "2017": "data/recoil_corrections/PuppiMETSys_2017.root",
                    "2018": "data/recoil_corrections/PuppiMETSys_2018.root",
                }
            ),
            "applyRecoilCorrections": SampleModifier({"wj": True}, default=False),
            "apply_recoil_resolution_systematic": False,
            "apply_recoil_response_systematic": False,
            "recoil_systematic_shift_up": False,
            "recoil_systematic_shift_down": False,
            "min_jetpt_met_propagation": 15,
        },
    )

    ## MT, MM channel trigger setup
    configuration.add_config_parameters(
        ["mt", "mm"],
        {
            "singlemoun_trigger": EraModifier(
                {
                    "2018": [
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
                }
            ),
            "cross_trigger": EraModifier(
                {
                    "2018": [
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
                }
            ),
        },
    )

    configuration.add_config_parameters(
        channels,
        {
            "ggHNNLOweightsRootfile": "data/htxs/NNLOPS_reweight.root",
            "ggH_generator": "powheg",
            "zptmass_file": EraModifier(
                {
                    "2016": "data/zpt/htt_scalefactors_legacy_2016.root",
                    "2017": "data/zpt/htt_scalefactors_legacy_2017.root",
                    "2018": "data/zpt/htt_scalefactors_legacy_2018.root",
                }
            ),
            "zptmass_functor": "zptmass_weight_nom",
            "zptmass_arguments": "z_gen_mass,z_gen_pt",
        },
    )
    configuration.add_producers(
        "global",
        [
            # RunLumiEventFilter,
            Lumi,
            MetFilter,
            PUweights,
            TauEnergyCorrection,
            GoodTaus,
            BaseMuons,
            BaseElectrons,
            JetEnergyCorrection,
            GoodJets,
            GoodBJets,
            DiLeptonVeto,
        ],
    )
    configuration.add_producers(
        "mm",
        [
            GoodMuons,
            MMPairSelection,
            GoodMMPairFilter,
            LVMu1,
            LVMu2,
            MMDiTauPairQuantities,
            JetCollection,
            BasicJetQuantities,
            BJetCollection,
            BasicBJetQuantities,
            MMGenDiTauPairQuantities,
            MuonIDIso_SF,
            MMGenerateSingleMuonTriggerFlags,
            LVMu1Uncorrected,
            LVMu2Uncorrected,
            MetCorrections,
            DiTauPairMETQuantities,
        ],
    )
    configuration.add_producers(
        "mt",
        [
            GoodMuons,
            MTPairSelection,
            GoodMTPairFilter,
            VetoMuons,
            ExtraMuonsVeto,
            ExtraElectronsVeto,
            LVMu1,
            LVTau2,
            MTDiTauPairQuantities,
            JetCollection,
            BasicJetQuantities,
            BJetCollection,
            BasicBJetQuantities,
            MTGenDiTauPairQuantities,
            MuonIDIso_SF,
            MTGenerateSingleMuonTriggerFlags,
            MTGenerateCrossTriggerFlags,
            LVMu1Uncorrected,
            LVTau2Uncorrected,
            MetCorrections,
            DiTauPairMETQuantities,
        ],
    )
    configuration.add_modification_rule(
        ["mt", "mm"], RemoveProducer(producers=[MuonIDIso_SF], samples=["data"])
    )
    configuration.add_modification_rule(
        "global",
        RemoveProducer(producers=[PUweights], samples=["data", "emb", "emb_mc"]),
    )
    configuration.add_modification_rule(
        ["mt", "mm"], AppendProducer(producers=[GGH_NNLO_Reweighting], samples=["ggh"])
    )
    configuration.add_modification_rule(
        ["mt", "mm"], AppendProducer(producers=[GGH_WG1_Uncertainties], samples=["ggh"])
    )
    configuration.add_modification_rule(
        ["mt", "mm"], AppendProducer(producers=[QQH_WG1_Uncertainties], samples=["qqh"])
    )
    configuration.add_modification_rule(
        ["mt", "mm"], AppendProducer(producers=[TopPtReweighting], samples=["ttbar"])
    )
    configuration.add_modification_rule(
        ["mt", "mm"], AppendProducer(producers=[ZPtMassReweighting], samples=["dy"])
    )
    # changes needed for data
    # global scope
    configuration.add_modification_rule(
        "global",
        AppendProducer(producers=[RenameJetsData], samples=["data", "emb", "emb_mc"]),
    )
    configuration.add_modification_rule(
        "global",
        AppendProducer(producers=[JSONFilter], samples=["data", "emb"]),
    )
    configuration.add_modification_rule(
        "global",
        RemoveProducer(
            producers=[JetEnergyCorrection], samples=["data", "emb", "emb_mc"]
        ),
    )
    # channel specific
    configuration.add_modification_rule(
        "mt",
        RemoveProducer(
            producers=[MTGenDiTauPairQuantities], samples=["data", "emb", "emb_mc"]
        ),
    )
    configuration.add_modification_rule(
        "mm",
        RemoveProducer(
            producers=[MMGenDiTauPairQuantities], samples=["data", "emb", "emb_mc"]
        ),
    )

    configuration.add_outputs(
        "mt",
        [
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
            q.pt_vis,
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
            q.met,
            q.metphi,
            q.metSumEt,
            q.metcov00,
            q.metcov01,
            q.metcov10,
            q.metcov11,
            MTGenerateSingleMuonTriggerFlags.output_group,
            MTGenerateCrossTriggerFlags.output_group,
            q.pzetamissvis,
            q.mTdileptonMET,
            q.mt_1,
            q.mt_2,
            q.pt_tt,
            q.pt_ttjj,
            q.mt_tot,
        ],
    )
    configuration.add_outputs(
        "mm",
        [
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
            q.pt_vis,
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
            q.idWeight_1,
            q.isoWeight_1,
            q.met,
            q.metphi,
            q.metSumEt,
            q.metcov00,
            q.metcov01,
            q.metcov10,
            q.metcov11,
            MMGenerateSingleMuonTriggerFlags.output_group,
            q.pzetamissvis,
            q.mTdileptonMET,
            q.mt_1,
            q.mt_2,
            q.pt_tt,
            q.pt_ttjj,
            q.mt_tot,
        ],
    )
    configuration.optimize()
    configuration.add_shift(
    SystematicShift(
    name="tauES_1prong0pizeroDown",
    shift_config={"global": {"tau_ES_shift_DM0": 0.998}},
    producers={"global": TauPtCorrection},
    ignore_producers={"mt": [LVMu1, VetoMuons]},
    ))
    # if "data" not in sample and "emb" not in sample:
    #     for scope in config["output"].keys():
    #         config["output"][scope].extend(
    #             [
    #                 nanoAOD.HTXS_Higgs_pt,
    #                 nanoAOD.HTXS_njets30,
    #                 nanoAOD.HTXS_stage_0,
    #                 nanoAOD.HTXS_stage1_2_cat_pTjet30GeV,
    #             ]
    #         )

    # for modifier in config["producer_modifiers"]:
    #     modifier.apply(sample, config["producers"], config["output"])
    # ResolveSampleDependencies(config, sample)
    # ResolveEraDependencies(config, era)
    # OptimizeProducerOrdering(config)
    # ValidateOutputs(config)
    # # remove channels that are not envoked
    # if "auto" not in channels:
    #     available_channels = [x for x in config["producers"] if x != "global"]
    #     for channel in available_channels:
    #         if channel not in channels:
    #             del config["producers"][channel]
    #             del config["output"][channel]

    # AddSystematicShift(
    #     config,
    #     "tauES_1prong0pizeroUp",
    #     {"global": {"tau_ES_shift_DM0": 1.002}},
    #     [[TauPtCorrection, "global"]],
    #     sanitize_producers=[[LVMu1, "mt"], [VetoMuons, "mt"]],
    # )
    # AddSystematicShift(
    #     config,
    #     "tauES_1prong0pizeroDown",
    #     {"global": {"tau_ES_shift_DM0": 0.998}},
    #     [[TauPtCorrection, "global"]],
    #     sanitize_producers=[[LVMu1, "mt"], [VetoMuons, "mt"]],
    # )
    # # Add MET shifts
    # if "data" not in sample and "emb" not in sample:
    #     SystematicShiftByInputQuantity(
    #         config,
    #         "metUnclusteredEnUp",
    #         {
    #             nanoAOD.MET_pt: "PuppiMET_ptUnclusteredUp",
    #             nanoAOD.MET_phi: "PuppiMET_phiUnclusteredUp",
    #         },
    #     )
    #     SystematicShiftByInputQuantity(
    #         config,
    #         "metUnclusteredEnDown",
    #         {
    #             nanoAOD.MET_pt: "PuppiMET_ptUnclusteredDown",
    #             nanoAOD.MET_phi: "PuppiMET_phiUnclusteredDown",
    #         },
    #     )
    #     # MET Recoil Shifts
    #     AddSystematicShift(
    #         config,
    #         "metRecoilResponseUp",
    #         {
    #             "mt": {
    #                 "apply_recoil_resolution_systematic": False,
    #                 "apply_recoil_response_systematic": True,
    #                 "recoil_systematic_shift_up": True,
    #                 "recoil_systematic_shift_down": False,
    #             }
    #         },
    #         [[ApplyRecoilCorrections, "mt"]],
    #     )
    #     AddSystematicShift(
    #         config,
    #         "metRecoilResponseDown",
    #         {
    #             "mt": {
    #                 "apply_recoil_resolution_systematic": False,
    #                 "apply_recoil_response_systematic": True,
    #                 "recoil_systematic_shift_up": False,
    #                 "recoil_systematic_shift_down": True,
    #             }
    #         },
    #         [[ApplyRecoilCorrections, "mt"]],
    #     )
    #     AddSystematicShift(
    #         config,
    #         "metRecoilResolutionUp",
    #         {
    #             "mt": {
    #                 "apply_recoil_resolution_systematic": True,
    #                 "apply_recoil_response_systematic": False,
    #                 "recoil_systematic_shift_up": True,
    #                 "recoil_systematic_shift_down": False,
    #             }
    #         },
    #         [[ApplyRecoilCorrections, "mt"]],
    #     )
    #     AddSystematicShift(
    #         config,
    #         "metRecoilResolutionDown",
    #         {
    #             "mt": {
    #                 "apply_recoil_resolution_systematic": True,
    #                 "apply_recoil_response_systematic": False,
    #                 "recoil_systematic_shift_up": False,
    #                 "recoil_systematic_shift_down": True,
    #             }
    #         },
    #         [[ApplyRecoilCorrections, "mt"]],
    #     )
    # # Jet energy resolution
    # shift_dict = {"JE_reso_shift": 1}
    # AddSystematicShift(
    #     config, "jerUncUp", {"global": shift_dict}, [[JetEnergyCorrection, "global"]]
    # )
    # shift_dict = {"JE_reso_shift": -1}
    # AddSystematicShift(
    #     config, "jerUncDown", {"global": shift_dict}, [[JetEnergyCorrection, "global"]]
    # )
    configuration.validate()
    configuration.report()
    return configuration.dump_dict()
