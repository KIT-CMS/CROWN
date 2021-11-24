from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List

import code_generation.producers.electrons as electrons
import code_generation.producers.event as event
import code_generation.producers.genparticles as genparticles
import code_generation.producers.jets as jets
import code_generation.producers.met as met
import code_generation.producers.muons as muons
import code_generation.producers.pairquantities as pairquantities
import code_generation.producers.pairselection as pairselection
import code_generation.producers.scalefactors as scalefactors
import code_generation.producers.taus as taus
import code_generation.quantities.nanoAOD as nanoAOD
import code_generation.quantities.output as q
from code_generation.configuration import Configuration
from code_generation.modifiers import EraModifier, SampleModifier
from code_generation.rules import AppendProducer, RemoveProducer
from code_generation.systematics import SystematicShift, SystematicShiftByQuantity


def build_config(
    era: str,
    sample: str,
    channels: List[str],
    shifts: List[str],
    available_sample_types: List[str],
    available_eras: List[str],
    available_channels: List[str],
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
    # MT/MM channel Muon selection
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
            event.Lumi,
            event.MetFilter,
            event.PUweights,
            taus.TauEnergyCorrection,
            taus.GoodTaus,
            muons.BaseMuons,
            electrons.BaseElectrons,
            jets.JetEnergyCorrection,
            jets.GoodJets,
            jets.GoodBJets,
            event.DiLeptonVeto,
        ],
    )
    configuration.add_producers(
        "mm",
        [
            muons.GoodMuons,
            pairselection.MMPairSelection,
            pairselection.GoodMMPairFilter,
            pairselection.LVMu1,
            pairselection.LVMu2,
            pairselection.LVMu1Uncorrected,
            pairselection.LVMu2Uncorrected,
            pairquantities.MMDiTauPairQuantities,
            jets.JetCollection,
            jets.BasicJetQuantities,
            jets.BJetCollection,
            jets.BasicBJetQuantities,
            genparticles.MMGenDiTauPairQuantities,
            scalefactors.MuonIDIso_SF,
            met.MetCorrections,
            pairquantities.DiTauPairMETQuantities,
        ],
    )
    configuration.add_producers(
        "mt",
        [
            muons.GoodMuons,
            pairselection.MTPairSelection,
            pairselection.GoodMTPairFilter,
            muons.VetoMuons,
            muons.ExtraMuonsVeto,
            electrons.ExtraElectronsVeto,
            pairselection.LVMu1,
            pairselection.LVTau2,
            pairselection.LVMu1Uncorrected,
            pairselection.LVTau2Uncorrected,
            pairquantities.MTDiTauPairQuantities,
            jets.JetCollection,
            jets.BasicJetQuantities,
            jets.BJetCollection,
            jets.BasicBJetQuantities,
            genparticles.MTGenDiTauPairQuantities,
            scalefactors.MuonIDIso_SF,
            met.MetCorrections,
            pairquantities.DiTauPairMETQuantities,
        ],
    )
    configuration.add_modification_rule(
        ["mt", "mm"],
        RemoveProducer(producers=scalefactors.MuonIDIso_SF, samples=["data"]),
    )
    configuration.add_modification_rule(
        "global",
        RemoveProducer(producers=event.PUweights, samples=["data", "emb", "emb_mc"]),
    )
    configuration.add_modification_rule(
        ["mt", "mm"],
        AppendProducer(
            producers=[event.GGH_NNLO_Reweighting, event.GGH_WG1_Uncertainties],
            samples=["ggh"],
        ),
    )
    configuration.add_modification_rule(
        ["mt", "mm"],
        AppendProducer(producers=event.QQH_WG1_Uncertainties, samples=["qqh"]),
    )
    configuration.add_modification_rule(
        ["mt", "mm"],
        AppendProducer(producers=event.TopPtReweighting, samples=["ttbar"]),
    )
    configuration.add_modification_rule(
        ["mt", "mm"],
        AppendProducer(producers=event.ZPtMassReweighting, samples=["dy"]),
    )
    # changes needed for data
    # global scope
    configuration.add_modification_rule(
        "global",
        AppendProducer(
            producers=jets.RenameJetsData, samples=["data", "emb", "emb_mc"]
        ),
    )
    configuration.add_modification_rule(
        "global",
        AppendProducer(producers=event.JSONFilter, samples=["data", "emb"]),
    )
    configuration.add_modification_rule(
        "global",
        RemoveProducer(
            producers=jets.JetEnergyCorrection, samples=["data", "emb", "emb_mc"]
        ),
    )
    # channel specific
    configuration.add_modification_rule(
        "mt",
        RemoveProducer(
            producers=[genparticles.MTGenDiTauPairQuantities],
            samples="data",
        ),
    )
    configuration.add_modification_rule(
        "mm",
        RemoveProducer(
            producers=[genparticles.MMGenDiTauPairQuantities],
            samples="data",
        ),
    )

    configuration.add_outputs(
        ["mt", "mm"],
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
            q.met_uncorrected,
            q.metphi_uncorrected,
            q.pfmet_uncorrected,
            q.pfmetphi_uncorrected,
            q.metSumEt,
            q.metcov00,
            q.metcov01,
            q.metcov10,
            q.metcov11,
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
        "mt",
        [
            q.taujet_pt_2,
            q.gen_taujet_pt_2,
            q.decaymode_2,
            q.gen_match_2,
            q.muon_veto_flag,
            q.dimuon_veto,
            q.electron_veto_flag,
        ],
    )
    #########################
    # Jet energy resolution
    #########################
    configuration.add_shift(
        SystematicShift(
            name="jerUncUp",
            shift_config={
                "global": {"JE_reso_shift": 1},
            },
            producers={"global": jets.JetEnergyCorrection},
        )
    )
    configuration.add_shift(
        SystematicShift(
            name="jerUncDown",
            shift_config={
                "global": {"JE_reso_shift": -1},
            },
            producers={"global": jets.JetEnergyCorrection},
        )
    )
    configuration.add_shift(
        SystematicShiftByQuantity(
            name="metUnclusteredEnUp",
            quantity_change={
                nanoAOD.MET_pt: "PuppiMET_ptUnclusteredUp",
                nanoAOD.MET_phi: "PuppiMET_phiUnclusteredUp",
            },
            scopes=["et", "mt", "tt", "em", "ee", "mm"],
        )
    )
    configuration.add_shift(
        SystematicShiftByQuantity(
            name="metUnclusteredEnDown",
            quantity_change={
                nanoAOD.MET_pt: "PuppiMET_ptUnclusteredDown",
                nanoAOD.MET_phi: "PuppiMET_phiUnclusteredDown",
            },
            scopes=["et", "mt", "tt", "em", "ee", "mm"],
        )
    )

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.dump_dict()
