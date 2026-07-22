from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, Filter
from code_generation.helpers import defaults

####################
# Set of producers used for contruction of MT good pairs and the coressponding lorentz vectors
####################

with defaults(scopes=["mt"]):
    MTPairSelection = Producer(
        call="ditau_pairselection::mutau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
        input=[
            q.Tau_pt_corrected,
            nanoAOD.Tau_eta,
            nanoAOD.Tau_phi,
            nanoAOD.Tau_mass,
            nanoAOD.Tau_rawDeepTau2017v2p1VSjet,
            nanoAOD.Muon_pt,
            nanoAOD.Muon_eta,
            nanoAOD.Muon_phi,
            nanoAOD.Muon_mass,
            nanoAOD.Muon_pfRelIso04_all,
            q.good_muons_mask,
            q.good_taus_mask,
        ],
        output=[q.dileptonpair],
    )

    GoodMTPairFlag = Producer(
        call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
        input=[q.dileptonpair],
        output=[],
    )

    GoodMTPairFilter = Filter(
        call='event::filter::Flags({df}, "GoodMuTauPairs", {input}, "any_of")',
        input=[],
        subproducers=[GoodMTPairFlag],
    )
with defaults(scopes=["mm"]):
    MMPairSelection = Producer(
        call="ditau_pairselection::mumu::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
        input=[
            nanoAOD.Muon_pt,
            nanoAOD.Muon_eta,
            nanoAOD.Muon_phi,
            nanoAOD.Muon_mass,
            q.good_muons_mask,
        ],
        output=[q.dileptonpair],
    )
    ZMMPairSelection = Producer(
        call="ditau_pairselection::mumu::ZBosonPairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
        input=[
            nanoAOD.Muon_pt,
            nanoAOD.Muon_eta,
            nanoAOD.Muon_phi,
            nanoAOD.Muon_mass,
            q.good_muons_mask,
        ],
        output=[q.dileptonpair],
    )

    GoodMMPairFlag = Producer(
        call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
        input=[q.dileptonpair],
        output=[],
    )

    GoodMMPairFilter = Filter(
        call='event::filter::Flags({df}, "GoodMuMuPairs", {input}, "any_of")',
        input=[],
        subproducers=[GoodMMPairFlag],
    )
with defaults(scopes=["et"]):
    ETPairSelection = Producer(
        call="ditau_pairselection::eltau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
        input=[
            q.Tau_pt_corrected,
            nanoAOD.Tau_eta,
            nanoAOD.Tau_phi,
            nanoAOD.Tau_mass,
            nanoAOD.Tau_rawDeepTau2017v2p1VSjet,
            nanoAOD.Electron_pt,
            nanoAOD.Electron_eta,
            nanoAOD.Electron_phi,
            nanoAOD.Electron_mass,
            nanoAOD.Electron_pfRelIso03_all,
            q.good_electrons_mask,
            q.good_taus_mask,
        ],
        output=[q.dileptonpair],
    )

    GoodETPairFlag = Producer(
        call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
        input=[q.dileptonpair],
        output=[],
    )

    GoodETPairFilter = Filter(
        call='event::filter::Flags({df}, "GoodElTauPairs", {input}, "any_of")',
        input=[],
        subproducers=[GoodETPairFlag],
    )
with defaults(scopes=["tt"]):
    TTPairSelection = Producer(
        call="ditau_pairselection::tautau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
        input=[
            q.Tau_pt_corrected,
            nanoAOD.Tau_eta,
            nanoAOD.Tau_phi,
            nanoAOD.Tau_mass,
            nanoAOD.Tau_rawDeepTau2017v2p1VSjet,
            q.good_taus_mask,
        ],
        output=[q.dileptonpair],
    )

    GoodTTPairFlag = Producer(
        call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
        input=[q.dileptonpair],
        output=[],
    )

    GoodTTPairFilter = Filter(
        call='event::filter::Flags({df}, "GoodTauTauPairs", {input}, "any_of")',
        input=[],
        subproducers=[GoodTTPairFlag],
    )
with defaults(scopes=["em"]):
    EMPairSelection = Producer(
        call="ditau_pairselection::elmu::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
        input=[
            nanoAOD.Electron_pt,
            nanoAOD.Electron_eta,
            nanoAOD.Electron_phi,
            nanoAOD.Electron_mass,
            nanoAOD.Electron_pfRelIso03_all,
            nanoAOD.Muon_pt,
            nanoAOD.Muon_eta,
            nanoAOD.Muon_phi,
            nanoAOD.Muon_mass,
            nanoAOD.Muon_pfRelIso04_all,
            q.good_electrons_mask,
            q.good_muons_mask,
        ],
        output=[q.dileptonpair],
    )

    GoodEMPairFlag = Producer(
        call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
        input=[q.dileptonpair],
        output=[],
    )

    GoodEMPairFilter = Filter(
        call='event::filter::Flags({df}, "GoodElMuPairs", {input}, "any_of")',
        input=[],
        subproducers=[GoodEMPairFlag],
    )

with defaults(call="lorentzvector::Build({df}, {output}, {input}, 0)"):
    with defaults(output=[q.p4_1]):
        LVMu1 = Producer(
            input=[
                nanoAOD.Muon_pt,
                nanoAOD.Muon_eta,
                nanoAOD.Muon_phi,
                nanoAOD.Muon_mass,
                q.dileptonpair,
            ],
            scopes=["mt", "mm"],
        )
        LVEl1 = Producer(
            input=[
                nanoAOD.Electron_pt,
                nanoAOD.Electron_eta,
                nanoAOD.Electron_phi,
                nanoAOD.Electron_mass,
                q.dileptonpair,
            ],
            scopes=["et", "ee", "em"],
        )
        LVTau1 = Producer(
            input=[
                q.Tau_pt_corrected,
                nanoAOD.Tau_eta,
                nanoAOD.Tau_phi,
                q.Tau_mass_corrected,
                q.dileptonpair,
            ],
            scopes=["tt"],
        )
    with defaults(output=[q.p4_1_uncorrected]):
        LVMu1Uncorrected = Producer(
            input=[
                nanoAOD.Muon_pt,
                nanoAOD.Muon_eta,
                nanoAOD.Muon_phi,
                nanoAOD.Muon_mass,
                q.dileptonpair,
            ],
            scopes=["mt", "mm"],
        )
        LVEl1Uncorrected = Producer(
            input=[
                nanoAOD.Electron_pt,
                nanoAOD.Electron_eta,
                nanoAOD.Electron_phi,
                nanoAOD.Electron_mass,
                q.dileptonpair,
            ],
            scopes=["em", "et", "ee"],
        )
        LVTau1Uncorrected = Producer(
            input=[
                nanoAOD.Tau_pt,
                nanoAOD.Tau_eta,
                nanoAOD.Tau_phi,
                nanoAOD.Tau_mass,
                q.dileptonpair,
            ],
            scopes=["tt"],
        )
with defaults(call="lorentzvector::Build({df}, {output}, {input}, 1)"):
    with defaults(output=[q.p4_2]):
        LVMu2 = Producer(
            input=[
                nanoAOD.Muon_pt,
                nanoAOD.Muon_eta,
                nanoAOD.Muon_phi,
                nanoAOD.Muon_mass,
                q.dileptonpair,
            ],
            scopes=["mm", "em"],
        )
        LVEl2 = Producer(
            input=[
                nanoAOD.Electron_pt,
                nanoAOD.Electron_eta,
                nanoAOD.Electron_phi,
                nanoAOD.Electron_mass,
                q.dileptonpair,
            ],
            scopes=["ee"],
        )
        LVTau2 = Producer(
            input=[
                q.Tau_pt_corrected,
                nanoAOD.Tau_eta,
                nanoAOD.Tau_phi,
                q.Tau_mass_corrected,
                q.dileptonpair,
            ],
            scopes=["mt", "et", "tt"],
        )
    with defaults(output=[q.p4_2_uncorrected]):
        LVMu2Uncorrected = Producer(
            input=[
                nanoAOD.Muon_pt,
                nanoAOD.Muon_eta,
                nanoAOD.Muon_phi,
                nanoAOD.Muon_mass,
                q.dileptonpair,
            ],
            scopes=["mm", "em"],
        )
        LVEl2Uncorrected = Producer(
            input=[
                nanoAOD.Electron_pt,
                nanoAOD.Electron_eta,
                nanoAOD.Electron_phi,
                nanoAOD.Electron_mass,
                q.dileptonpair,
            ],
            scopes=["ee"],
        )
        LVTau2Uncorrected = Producer(
            input=[
                nanoAOD.Tau_pt,
                nanoAOD.Tau_eta,
                nanoAOD.Tau_phi,
                nanoAOD.Tau_mass,
                q.dileptonpair,
            ],
            scopes=["mt", "et", "tt"],
        )
