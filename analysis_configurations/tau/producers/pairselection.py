import config.tau.quantities.output as q
import config.tau.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, Filter

####################
# Set of producers used for contruction of MT good pairs and the coressponding lorentz vectors
####################

MTPairSelection = Producer(
    name="MTPairSelection",
    call="pairselection::mutau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
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
    output=[q.ditaupair],
    scopes=["mt"],
)

GoodMTPairFlag = Producer(
    name="GoodMTPairFlag",
    call="pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.ditaupair],
    output=[],
    scopes=["mt"],
)

GoodMTPairFilter = Filter(
    name="GoodMTPairFilter",
    call='basefunctions::FilterFlagsAny({df}, "GoodMuTauPairs", {input})',
    input=[],
    scopes=["mt"],
    subproducers=[GoodMTPairFlag],
)

MMPairSelection = Producer(
    name="MMPairSelection",
    call="pairselection::mumu::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        q.good_muons_mask,
    ],
    output=[q.ditaupair],
    scopes=["mm"],
)
ZMMPairSelection = Producer(
    name="MMPairSelection",
    call="pairselection::mumu::ZBosonPairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        q.good_muons_mask,
    ],
    output=[q.ditaupair],
    scopes=["mm"],
)

GoodMMPairFlag = Producer(
    name="GoodMMPairFlag",
    call="pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.ditaupair],
    output=[],
    scopes=["mm"],
)

GoodMMPairFilter = Filter(
    name="GoodMMPairFilter",
    call='basefunctions::FilterFlagsAny({df}, "GoodMuMuPairs", {input})',
    input=[],
    scopes=["mm"],
    subproducers=[GoodMMPairFlag],
)

ETPairSelection = Producer(
    name="ETPairSelection",
    call="pairselection::eltau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
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
    output=[q.ditaupair],
    scopes=["et"],
)

GoodETPairFlag = Producer(
    name="GoodETPairFlag",
    call="pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.ditaupair],
    output=[],
    scopes=["et"],
)

GoodETPairFilter = Filter(
    name="GoodETPairFilter",
    call='basefunctions::FilterFlagsAny({df}, "GoodElTauPairs", {input})',
    input=[],
    scopes=["et"],
    subproducers=[GoodETPairFlag],
)

####################
## TauTau Pair Selection
####################
TTPairSelection = Producer(
    name="TTPairSelection",
    call="pairselection::tautau::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
    input=[
        q.Tau_pt_corrected,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        nanoAOD.Tau_mass,
        nanoAOD.Tau_IDraw,
        q.good_taus_mask,
    ],
    output=[q.ditaupair],
    scopes=["tt"],
)

GoodTTPairFlag = Producer(
    name="GoodTTPairFlag",
    call="pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.ditaupair],
    output=[],
    scopes=["tt"],
)

GoodTTPairFilter = Filter(
    name="GoodTTPairFilter",
    call='basefunctions::FilterFlagsAny({df}, "GoodTauTauPairs", {input})',
    input=[],
    scopes=["tt"],
    subproducers=[GoodTTPairFlag],
)
####################
## ElMu Pair Selection
####################

EMPairSelection = Producer(
    name="EMPairSelection",
    call="pairselection::elmu::PairSelection({df}, {input_vec}, {output}, {pairselection_min_dR})",
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
    output=[q.ditaupair],
    scopes=["em"],
)

GoodEMPairFlag = Producer(
    name="GoodEMPairFlag",
    call="pairselection::flagGoodPairs({df}, {output}, {input})",
    input=[q.ditaupair],
    output=[],
    scopes=["em"],
)

GoodEMPairFilter = Filter(
    name="GoodEMPairFilter",
    call='basefunctions::FilterFlagsAny({df}, "GoodElMuPairs", {input})',
    input=[],
    scopes=["em"],
    subproducers=[GoodEMPairFlag],
)


LVMu1 = Producer(
    name="LVMu1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
    ],
    output=[q.p4_1],
    scopes=["mt", "mm"],
)
LVMu2 = Producer(
    name="LVMu2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
    ],
    output=[q.p4_2],
    scopes=["mm", "em"],
)
LVEl1 = Producer(
    name="LVEl1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
    ],
    output=[q.p4_1],
    scopes=["et", "ee", "em"],
)
LVEl2 = Producer(
    name="LVEl2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
    ],
    output=[q.p4_2],
    scopes=["ee"],
)
LVTau1 = Producer(
    name="LVTau1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.ditaupair,
        q.Tau_pt_corrected,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        q.Tau_mass_corrected,
    ],
    output=[q.p4_1],
    scopes=["tt"],
)
LVTau2 = Producer(
    name="LVTau2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.ditaupair,
        q.Tau_pt_corrected,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        q.Tau_mass_corrected,
    ],
    output=[q.p4_2],
    scopes=["mt", "et", "tt"],
)
## uncorrected versions of all particles, used for MET propagation
LVMu1Uncorrected = Producer(
    name="LVMu1Uncorrected",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
    ],
    output=[q.p4_1_uncorrected],
    scopes=["mt", "mm"],
)
LVMu2Uncorrected = Producer(
    name="LVMu2Uncorrected",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
    ],
    output=[q.p4_2_uncorrected],
    scopes=["mm", "em"],
)
LVEl1Uncorrected = Producer(
    name="LVEl1Uncorrected",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
    ],
    output=[q.p4_1_uncorrected],
    scopes=["em", "et", "ee"],
)
LVEl2Uncorrected = Producer(
    name="LVEl2Uncorrected",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Electron_pt,
        nanoAOD.Electron_eta,
        nanoAOD.Electron_phi,
        nanoAOD.Electron_mass,
    ],
    output=[q.p4_2_uncorrected],
    scopes=["ee"],
)
LVTau1Uncorrected = Producer(
    name="LVTau1Uncorrected",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Tau_pt,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        nanoAOD.Tau_mass,
    ],
    output=[q.p4_1_uncorrected],
    scopes=["tt"],
)
LVTau2Uncorrected = Producer(
    name="LVTau2Uncorrected",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.ditaupair,
        nanoAOD.Tau_pt,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_phi,
        nanoAOD.Tau_mass,
    ],
    output=[q.p4_2_uncorrected],
    scopes=["mt", "et", "tt"],
)
