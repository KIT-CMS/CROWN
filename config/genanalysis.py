from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List
from code_generation.modifiers import SampleModifier

import code_generation.producers.event as event
import code_generation.producers.met as met
import code_generation.producers.muons as muons
import code_generation.producers.pairquantities as pairquantities
import code_generation.producers.genparticles as genparticles
import code_generation.producers.pairselection as pairselection
import code_generation.producers.embedding as emb
import code_generation.quantities.nanoAOD as nanoAOD
import code_generation.quantities.output as q
from code_generation.configuration import Configuration
from code_generation.rules import AppendProducer


def build_config(
    era: str,
    sample: str,
    channels: List[str],
    shifts: List[str],
    available_sample_types: List[str],
    available_eras: List[str],
    available_channels: List[str],
):

    if sample == "data":
        print("WARNING: no genparticles available in data, what are you doing ??")
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
            "RunLumiEventFilter_Quantities": ["event", "luminosityBlock"],
            "RunLumiEventFilter_Quantity_Types": ["ULong64_t", "UInt_t"],
            "RunLumiEventFilter_Selections": ["3", "318"],
            "min_muon_pt": 10.0,
            "max_muon_eta": 2.4,
            "max_muon_dxy": 0.5,
            "max_muon_dz": 0.5,
            "muon_id": "Muon_looseId",
            "muon_iso_cut": 2.5,
        },
    )
    ###### Channel Specifics ######
    # MT/MM channel Muon selection
    configuration.add_config_parameters(
        ["mm"],
        {
            "muon_index_in_pair": 0,
            "second_muon_index_in_pair": 1,
            "muon_iso_cut": 0.15,
            "min_muon_pt": 23.0,
            "max_muon_eta": 2.4,
            "max_muon_dxy": 0.045,
            "max_muon_dz": 0.2,
            "muon_id": "Muon_mediumId",
            "truegen_mother_pdgid": SampleModifier(
                {"emb_mc": 23, "tt": 6, "vv": 24}, default=None
            ),
            "truegen_daughter_1_pdgid": 13,
            "truegen_daugher_2_pdgid": 13,
        },
    )

    configuration.add_producers(
        "global",
        [
            # event.RunLumiEventFilter,
            event.Lumi,
            muons.BaseMuons,
        ],
    )
    configuration.add_producers(
        "mm",
        [
            genparticles.MMTrueGenDiTauPairQuantities,
            met.UncorrectedMet,
            muons.GoodMuons,
            muons.VetoMuons,
            muons.VetoSecondMuon,
            muons.ExtraMuonsVeto,
            pairselection.ZMMPairSelection,
            pairselection.GoodMMPairFilter,
            pairselection.LVMu1,
            pairselection.LVMu2,
            pairselection.LVMu1Uncorrected,
            pairselection.LVMu2Uncorrected,
            pairquantities.MMDiTauPairQuantities,
        ],
    )

    configuration.add_outputs(
        ["mm"],
        [
            nanoAOD.run,
            q.lumi,
            nanoAOD.event,
            q.gen_pt_1,
            q.gen_eta_1,
            q.gen_phi_1,
            q.gen_mass_1,
            q.gen_pt_2,
            q.gen_eta_2,
            q.gen_phi_2,
            q.gen_mass_2,
            q.gen_m_vis,
            q.pt_1,
            q.pt_2,
            q.eta_1,
            q.eta_2,
            q.phi_1,
            q.phi_2,
            q.m_vis,
            q.muon_veto_flag,
        ],
    )

    configuration.add_modification_rule(
        "global",
        AppendProducer(producers=emb.EmbeddingQuantities, samples=["emb", "emb_mc"]),
    )

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.dump_dict()
