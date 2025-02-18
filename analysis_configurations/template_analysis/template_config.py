from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List

from .producers import event as event
from .producers import genparticles as genparticles
from .producers import muons as muons
from .producers import pairquantities as pairquantities
from .producers import pairselection as pairselection
from .producers import scalefactors as scalefactors
from .quantities import nanoAOD as nanoAOD
from .quantities import output as q
from code_generation.configuration import Configuration
from code_generation.modifiers import EraModifier
from code_generation.rules import RemoveProducer, AppendProducer
from code_generation.systematics import SystematicShift


def build_config(
    era: str,
    sample: str,
    scopes: List[str],
    shifts: List[str],
    available_sample_types: List[str],
    available_eras: List[str],
    available_scopes: List[str],
):
    configuration = Configuration(
        era,
        sample,
        scopes,
        shifts,
        available_sample_types,
        available_eras,
        available_scopes,
    )

    configuration.add_config_parameters(
        "global",
        {
            "golden_json_file": EraModifier(
                {
                    "2016preVFP": "data/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
                    "2016postVFP": "data/golden_json/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
                    "2017": "data/golden_json/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
                    "2018": "data/golden_json/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
                }
            ),
            "met_filters": EraModifier(
                {
                    "2016preVFP": [
                        "Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        "Flag_BadPFMuonDzFilter",
                        "Flag_eeBadScFilter",
                    ],
                    "2016postVFP": [
                        "Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        "Flag_BadPFMuonDzFilter",
                        "Flag_eeBadScFilter",
                    ],
                    "2017": [
                        "Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        "Flag_BadPFMuonDzFilter",
                        "Flag_eeBadScFilter",
                        "Flag_ecalBadCalibFilter",
                    ],
                    "2018": [
                        "Flag_goodVertices",
                        "Flag_globalSuperTightHalo2016Filter",
                        "Flag_HBHENoiseFilter",
                        "Flag_HBHENoiseIsoFilter",
                        "Flag_EcalDeadCellTriggerPrimitiveFilter",
                        "Flag_BadPFMuonFilter",
                        "Flag_BadPFMuonDzFilter",
                        "Flag_eeBadScFilter",
                        "Flag_ecalBadCalibFilter",
                    ],
                }
            ),
        },
    )
    # muon base selection applied before given to a specific scope
    configuration.add_config_parameters(
        "global",
        {
            "muon_min_pt": 10.0,
            "muon_max_eta": 2.4,
            "muon_max_dxy": 0.045,
            "muon_max_dz": 0.2,
            "muon_id": "Muon_mediumId",
            "muon_max_iso": 0.3,
        },
    )
    # MM scope Muon selection: additional cut on top of base selection in "global" scope
    configuration.add_config_parameters(
        ["mm"],
        {
            
            "muon_first_index_in_pair": 0,
            "muon_second_index_in_pair": 1,
            "muon_min_pt": 25.0,
            "muon_max_eta": 2.1,
            "muon_max_iso": 0.15,
        },
    )
    # Muon scale factors configuration based on jsonpog correctionlib files
    configuration.add_config_parameters(
        ["mm"],
        {
            "muon_sf_file": EraModifier(
                {
                    "2016preVFP": "data/jsonpog-integration/POG/MUO/2016preVFP_UL/muon_Z.json.gz",
                    "2016postVFP": "data/jsonpog-integration/POG/MUO/2016postVFP_UL/muon_Z.json.gz",
                    "2017": "data/jsonpog-integration/POG/MUO/2017_UL/muon_Z.json.gz",
                    "2018": "data/jsonpog-integration/POG/MUO/2018_UL/muon_Z.json.gz",
                }
            ),
            # relevant inputs for scale factor evaluation
            "muon_reco_sf_name": "NUM_TrackerMuons_DEN_genTracks",
            "muon_id_sf_name": "NUM_MediumID_DEN_TrackerMuons",
            "muon_iso_sf_name": "NUM_LooseRelIso_DEN_MediumID",
            "muon_sf_varation": "nominal",  # "systup"/"systdown" are up/down variations
        },
    )
    ## all scopes settings: applied in all chosen final scopes the same way
    configuration.add_config_parameters(
        scopes,
        {
            "mm_pair_min_deltaR": 0.5,
        },
    )

    configuration.add_producers(
        "global",
        [
            event.SampleFlags,
            event.Lumi,
            event.MetFilter,
            muons.BaseMuons,
        ],
    )
    configuration.add_producers(
        "mm",
        [
            muons.GoodMuons,
            muons.NumberOfGoodMuons,
            pairselection.ZMMPairSelection,
            pairselection.GoodMMPairFilter,
            pairselection.LVMu1,
            pairselection.LVMu2,
            pairquantities.MMDiTauPairQuantities,
            genparticles.MMGenDiTauPairQuantities,
            scalefactors.Muon_SFs,
        ],
    )

    configuration.add_outputs(
        "mm",
        [
            q.is_data,
            q.is_embedding,
            q.is_ttbar,
            q.is_dyjets,
            q.is_wjets,
            q.is_diboson,
            nanoAOD.run,
            q.lumi,
            nanoAOD.event,
            q.nmuons,
            q.pt_1,
            q.pt_2,
            q.eta_1,
            q.eta_2,
            q.phi_1,
            q.phi_2,
            q.mass_1,
            q.mass_2,
            q.q_1,
            q.q_2,
            q.mm_pair_mass,
            q.mm_pair_pt,
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
            q.gen_mm_pair_mass,
            q.reco_wgt_mu_1,
            q.reco_wgt_mu_2,
            q.id_wgt_mu_1,
            q.id_wgt_mu_2,
            q.iso_wgt_mu_1,
            q.iso_wgt_mu_2,
        ],
    )

    configuration.add_modification_rule(
        "mm",
        RemoveProducer(
            producers=[
                genparticles.MMGenDiTauPairQuantities,
                scalefactors.Muon_SFs,
            ],
            samples=["data"],
        ),
    )

    configuration.add_modification_rule(
        "global",
        AppendProducer(producers=event.JSONFilter, samples=["data", "embedding"]),
    )

    configuration.add_shift(
        SystematicShift(
            name="MuonIDUp",
            shift_config={"mm": {"muon_sf_varation": "systup"}},
            producers={
                "mm": [
                    scalefactors.Muon_1_ID_SF,
                    scalefactors.Muon_2_ID_SF,
                ]
            },
        )
    )
    configuration.add_shift(
        SystematicShift(
            name="MuonIDDown",
            shift_config={"mm": {"muon_sf_varation": "systdown"}},
            producers={
                "mm": [
                    scalefactors.Muon_1_ID_SF,
                    scalefactors.Muon_2_ID_SF,
                ]
            },
        )
    )

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()
