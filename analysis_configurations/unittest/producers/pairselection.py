from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, Filter

####################
# Set of producers used for contruction of MT good pairs and the coressponding lorentz vectors
####################

MTPairSelection = Producer(
    name="MTPairSelection",
    call="ditau_pairselection::mutau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        q.Tau_pt_corrected,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        nanoAOD.Tau_mass,
        nanoAOD.Tau_IDraw,
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        nanoAOD.Muon_iso,
        q.good_muons_mask,
        q.good_taus_mask,
    ],
    output=[q.dileptonpair],
    scopes=["mt"],
)

GoodMTPairFlag = Producer(
    name="GoodMTPairFlag",
    call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.dileptonpair],
    output=[],
    scopes=["mt"],
)

GoodMTPairFilter = Filter(
    name="GoodMTPairFilter",
    call='event::filter::Flags({df}, "GoodMuTauPairs", {input}, "any_of")',
    input=[],
    scopes=["mt"],
    subproducers=[GoodMTPairFlag],
)

MMPairSelection = Producer(
    name="MMPairSelection",
    call="ditau_pairselection::mumu::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        q.good_muons_mask,
    ],
    output=[q.dileptonpair],
    scopes=["mm"],
)
ZMMPairSelection = Producer(
    name="MMPairSelection",
    call="ditau_pairselection::mumu::ZBosonPairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        q.good_muons_mask,
    ],
    output=[q.dileptonpair],
    scopes=["mm"],
)

GoodMMPairFlag = Producer(
    name="GoodMMPairFlag",
    call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.dileptonpair],
    output=[],
    scopes=["mm"],
)

GoodMMPairFilter = Filter(
    name="GoodMMPairFilter",
    call='event::filter::Flags({df}, "GoodMuMuPairs", {input}, "any_of")',
    input=[],
    scopes=["mm"],
    subproducers=[GoodMMPairFlag],
)

ETPairSelection = Producer(
    name="ETPairSelection",
    call="ditau_pairselection::eltau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        q.Tau_pt_corrected,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        nanoAOD.Tau_mass,
        nanoAOD.Tau_IDraw,
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
        nanoAOD.Electron_iso,
        q.good_electrons_mask,
        q.good_taus_mask,
    ],
    output=[q.dileptonpair],
    scopes=["et"],
)

GoodETPairFlag = Producer(
    name="GoodETPairFlag",
    call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.dileptonpair],
    output=[],
    scopes=["et"],
)

GoodETPairFilter = Filter(
    name="GoodETPairFilter",
    call='event::filter::Flags({df}, "GoodElTauPairs", {input}, "any_of")',
    input=[],
    scopes=["et"],
    subproducers=[GoodETPairFlag],
)

####################
## TauTau Pair Selection
####################
TTPairSelection = Producer(
    name="TTPairSelection",
    call="ditau_pairselection::tautau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        q.Tau_pt_corrected,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        nanoAOD.Tau_mass,
        nanoAOD.Tau_IDraw,
        q.good_taus_mask,
    ],
    output=[q.dileptonpair],
    scopes=["tt"],
)

GoodTTPairFlag = Producer(
    name="GoodTTPairFlag",
    call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.dileptonpair],
    output=[],
    scopes=["tt"],
)

GoodTTPairFilter = Filter(
    name="GoodTTPairFilter",
    call='event::filter::Flags({df}, "GoodTauTauPairs", {input}, "any_of")',
    input=[],
    scopes=["tt"],
    subproducers=[GoodTTPairFlag],
)
####################
## ElMu Pair Selection
####################

EMPairSelection = Producer(
    name="EMPairSelection",
    call="ditau_pairselection::elmu::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
        nanoAOD.Electron_iso,
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        nanoAOD.Muon_iso,
        q.good_electrons_mask,
        q.good_muons_mask,
    ],
    output=[q.dileptonpair],
    scopes=["em"],
)

GoodEMPairFlag = Producer(
    name="GoodEMPairFlag",
    call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.dileptonpair],
    output=[],
    scopes=["em"],
)

GoodEMPairFilter = Filter(
    name="GoodEMPairFilter",
    call='event::filter::Flags({df}, "GoodElMuPairs", {input}, "any_of")',
    input=[],
    scopes=["em"],
    subproducers=[GoodEMPairFlag],
)


LVMu1 = Producer(
    name="LVMu1",
    call="lorentzvector::Build({df}, {output}, {input}, 0)",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        q.dileptonpair,
    ],
    output=[q.p4_1],
    scopes=["mt", "mm"],
)
LVMu2 = Producer(
    name="LVMu2",
    call="lorentzvector::Build({df}, {output}, {input}, 1)",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        q.dileptonpair,
    ],
    output=[q.p4_2],
    scopes=["mm", "em"],
)
LVEl1 = Producer(
    name="LVEl1",
    call="lorentzvector::Build({df}, {output}, {input}, 0)",
    input=[
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
        q.dileptonpair,
    ],
    output=[q.p4_1],
    scopes=["et", "ee", "em"],
)
LVEl2 = Producer(
    name="LVEl2",
    call="lorentzvector::Build({df}, {output}, {input}, 1)",
    input=[
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
        q.dileptonpair,
    ],
    output=[q.p4_2],
    scopes=["ee"],
)
LVTau1 = Producer(
    name="LVTau1",
    call="lorentzvector::Build({df}, {output}, {input}, 0)",
    input=[
        q.Tau_pt_corrected,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        q.Tau_mass_corrected,
        q.dileptonpair,
    ],
    output=[q.p4_1],
    scopes=["tt"],
)
LVTau2 = Producer(
    name="LVTau2",
    call="lorentzvector::Build({df}, {output}, {input}, 1)",
    input=[
        q.Tau_pt_corrected,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        q.Tau_mass_corrected,
        q.dileptonpair,
    ],
    output=[q.p4_2],
    scopes=["mt", "et", "tt"],
)
## uncorrected versions of all particles, used for MET propagation
LVMu1Uncorrected = Producer(
    name="LVMu1Uncorrected",
    call="lorentzvector::Build({df}, {output}, {input}, 0)",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        q.dileptonpair,
    ],
    output=[q.p4_1_uncorrected],
    scopes=["mt", "mm"],
)
LVMu2Uncorrected = Producer(
    name="LVMu2Uncorrected",
    call="lorentzvector::Build({df}, {output}, {input}, 1)",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        q.dileptonpair,
    ],
    output=[q.p4_2_uncorrected],
    scopes=["mm", "em"],
)
LVEl1Uncorrected = Producer(
    name="LVEl1Uncorrected",
    call="lorentzvector::Build({df}, {output}, {input}, 0)",
    input=[
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
        q.dileptonpair,
    ],
    output=[q.p4_1_uncorrected],
    scopes=["em", "et", "ee"],
)
LVEl2Uncorrected = Producer(
    name="LVEl2Uncorrected",
    call="lorentzvector::Build({df}, {output}, {input}, 1)",
    input=[
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
        q.dileptonpair,
    ],
    output=[q.p4_2_uncorrected],
    scopes=["ee"],
)
LVTau1Uncorrected = Producer(
    name="LVTau1Uncorrected",
    call="lorentzvector::Build({df}, {output}, {input}, 0)",
    input=[
        nanoAOD.Tau_pt,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        nanoAOD.Tau_mass,
        q.dileptonpair,
    ],
    output=[q.p4_1_uncorrected],
    scopes=["tt"],
)
LVTau2Uncorrected = Producer(
    name="LVTau2Uncorrected",
    call="lorentzvector::Build({df}, {output}, {input}, 1)",
    input=[
        nanoAOD.Tau_pt,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        nanoAOD.Tau_mass,
        q.dileptonpair,
    ],
    output=[q.p4_2_uncorrected],
    scopes=["mt", "et", "tt"],
)
