from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List
from code_generation.configuration import Configuration
from code_generation.systematics import SystematicShift
from .producers import jets as jets
from .producers import scalefactors as scalefactors


def add_jetVariations(
    configuration: Configuration, available_sample_types: List[str], era: str
):
    #########################
    # Jet energy resolution
    #########################
    configuration.add_shift(
        SystematicShift(
            name="jerUncUp",
            shift_config={
                "global": {"jet_jer_shift": '"up"'},
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jerUncDown",
            shift_config={
                "global": {"jet_jer_shift": '"down"'},
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    #########################
    # Jet energy scale - Total
    #########################
    JEC_sources = '{"Total"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncTotalUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                },
                ("mt", "et", "tt"): {"btag_sf_variation": "up_jes"},
            },
            producers={
                "global": {
                    jets.JetEnergyCorrection,
                },
                ("mt", "et", "tt"): {scalefactors.btagging_SF},
            },
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncTotalDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                },
                ("mt", "et", "tt"): {"btag_sf_variation": "up_jes"},
            },
            producers={
                "global": {
                    jets.JetEnergyCorrection,
                },
                ("mt", "et", "tt"): {scalefactors.btagging_SF},
            },
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    #########################
    # HEM 15/16 issue
    #########################
    if era == "2018":
        JEC_sources = '{"HEMIssue"}'
        configuration.add_shift(
            SystematicShift(
                name="jesUncHEMIssueUp",
                shift_config={
                    "global": {
                        "jet_jes_shift": 1,
                        "jet_jes_sources": JEC_sources,
                    }
                },
                producers={"global": jets.JetEnergyCorrection},
            ),
            samples=[
                sample
                for sample in available_sample_types
                if sample not in ["data", "embedding", "embedding_mc"]
            ],
        )
        configuration.add_shift(
            SystematicShift(
                name="jesUncHEMIssueDown",
                shift_config={
                    "global": {
                        "jet_jes_shift": -1,
                        "jet_jes_sources": JEC_sources,
                    }
                },
                producers={"global": jets.JetEnergyCorrection},
            ),
            samples=[
                sample
                for sample in available_sample_types
                if sample not in ["data", "embedding", "embedding_mc"]
            ],
        )
    #########################
    # Jet energy scale - individual
    #########################
    # JEC_sources = '{"AbsoluteStat"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncAbsoluteStatUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesAbsoluteStat",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncAbsoluteStatDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesAbsoluteStat",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"AbsoluteScale"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncAbsoluteScaleUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesAbsoluteScale",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncAbsoluteScaleDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesAbsoluteScale",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"AbsoluteMPFBias"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncAbsoluteMPFBiasUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesAbsoluteMPFBias",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncAbsoluteMPFBiasDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesAbsoluteMPFBias",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"Fragmentation"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncFragmentationUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesFragmentation",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncFragmentationDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesFragmentation",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"SinglePionECAL"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncSinglePionECALUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesSinglePionECAL",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncSinglePionECALDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesSinglePionECAL",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"SinglePionHCAL"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncSinglePionHCALUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesSinglePionHCAL",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncSinglePionHCALDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesSinglePionHCAL",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"FlavorQCD"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncFlavorQCDUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesFlavorQCD",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncFlavorQCDDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesFlavorQCD",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"TimePtEta"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncTimePtEtaUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesTimePtEta",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncTimePtEtaDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesTimePtEta",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeJEREC1"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeJEREC1Up",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeJEREC1",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeJEREC1Down",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeJEREC1",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeJEREC2"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeJEREC2Up",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeJEREC2",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeJEREC2Down",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeJEREC2",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeJERHF"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeJERHFUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeJERHF",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeJERHFDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeJERHF",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativePtBB"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativePtBBUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativePtBB",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativePtBBDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativePtBB",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativePtEC1"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativePtEC1Up",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativePtEC1",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativePtEC1Down",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativePtEC1",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativePtEC2"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativePtEC2Up",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativePtEC2",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativePtEC2Down",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativePtEC2",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativePtHF"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativePtHFUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativePtHF",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativePtHFDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativePtHF",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeBal"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeBalUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeBal",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeBalDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeBal",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeSample"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeSampleUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeSample",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeSampleDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeSample",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeFSR"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeFSRUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeFSR",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeFSRDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeFSR",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeStatFSR"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeStatFSRUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeStatFSR",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeStatFSRDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeStatFSR",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeStatEC"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeStatECUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeStatEC",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeStatECDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeStatEC",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"RelativeStatHF"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeStatHFUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesRelativeStatHF",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncRelativeStatHFDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesRelativeStatHF",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"PileUpDataMC"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpDataMCUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesPileUpDataMC",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpDataMCDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesPileUpDataMC",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"PileUpPtRef"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtRefUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesPileUpPtRef",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtRefDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesPileUpPtRef",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"PileUpPtBB"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtBBUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesPileUpPtBB",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtBBDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesPileUpPtBB",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"PileUpPtEC1"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtEC1Up",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesPileUpPtEC1",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtEC1Down",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesPileUpPtEC1",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"PileUpPtEC2"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtEC2Up",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesPileUpPtEC2",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtEC2Down",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesPileUpPtEC2",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    # JEC_sources = '{"PileUpPtHF"}'
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtHFUp",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": 1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "up_jesPileUpPtHF",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )
    # configuration.add_shift(
    #     SystematicShift(
    #         name="jesUncPileUpPtHFDown",
    #         shift_config={
    #             "global": {
    #                 "jet_jes_shift": -1,
    #                 "jet_jes_sources": JEC_sources,
    #                 "btag_sf_variation": "down_jesPileUpPtHF",
    #             }
    #         },
    #         producers={"global": jets.JetEnergyCorrection},
    #     ),
    #     samples=[
    #         sample
    #         for sample in available_sample_types
    #         if sample not in ["data", "embedding", "embedding_mc"]
    #     ],
    # )

    #########################
    # Jet energy scale - reduced set
    #########################
    JEC_sources = '{"SinglePionECAL", "SinglePionHCAL", "AbsoluteMPFBias", "AbsoluteScale", "Fragmentation", "PileUpDataMC", "RelativeFSR", "PileUpPtRef"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncAbsoluteUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncAbsoluteDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"AbsoluteStat", "TimePtEta", "RelativeStatFSR"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncAbsoluteYearUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncAbsoluteYearDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"FlavorQCD"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncFlavorQCDUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncFlavorQCDDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"PileUpPtEC1", "PileUpPtBB", "RelativePtBB"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncBBEC1Up",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncBBEC1Down",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"RelativeJEREC1", "RelativePtEC1", "RelativeStatEC"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncBBEC1YearUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncBBEC1YearDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"RelativePtHF", "PileUpPtHF", "RelativeJERHF"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncHFUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncHFDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"RelativeStatHF"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncHFYearUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncHFYearDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"PileUpPtEC2"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncEC2Up",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncEC2Down",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"RelativeJEREC2", "RelativePtEC2"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncEC2YearUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncEC2YearDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"RelativeBal"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncRelativeBalUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncRelativeBalDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    JEC_sources = '{"RelativeSample"}'
    configuration.add_shift(
        SystematicShift(
            name="jesUncRelativeSampleYearUp",
            shift_config={
                "global": {
                    "jet_jes_shift": 1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )
    configuration.add_shift(
        SystematicShift(
            name="jesUncRelativeSampleYearDown",
            shift_config={
                "global": {
                    "jet_jes_shift": -1,
                    "jet_jes_sources": JEC_sources,
                }
            },
            producers={"global": jets.JetEnergyCorrection},
        ),
        samples=[
            sample
            for sample in available_sample_types
            if sample not in ["data", "embedding", "embedding_mc"]
        ],
    )

    return configuration
