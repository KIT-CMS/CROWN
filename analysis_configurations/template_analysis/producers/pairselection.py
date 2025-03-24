from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, Filter

####################
# Set of producers used for contruction of MT good pairs and the coressponding lorentz vectors
####################

# mm pair with highest pt
MMPairSelection = Producer(
    name="MMPairSelection",
    call="ditau_pairselection::mumu::PairSelection({df}, {input_vec}, {output}, {mm_pair_min_deltaR})",
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
# mm pair closest to Z boson mass
ZMMPairSelection = Producer(
    name="MMPairSelection",
    call="ditau_pairselection::mumu::ZBosonPairSelection({df}, {input_vec}, {output}, {mm_pair_min_deltaR})",
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
    call='event::FilterFlagsAny({df}, "GoodMuMuPairs", {input})',
    input=[],
    scopes=["mm"],
    subproducers=[GoodMMPairFlag],
)


LVMu1 = Producer(
    name="LVMu1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.dileptonpair,
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
    ],
    output=[q.p4_1],
    scopes=["mm"],
)
LVMu2 = Producer(
    name="LVMu2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.dileptonpair,
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
    ],
    output=[q.p4_2],
    scopes=["mm"],
)

# Lorentz vectors for friend trees
LVMu1_friend = Producer(
    name="LVMu1_friend",
    call="lorentzvectors::build({df}, {input_vec}, {output})",
    input=[
        q.pt_1,
        q.eta_1,
        q.phi_1,
        q.mass_1,
    ],
    output=[q.p4_1],
    scopes=["mm"],
)
LVMu2_friend = Producer(
    name="LVMu2_friend",
    call="lorentzvectors::build({df}, {input_vec}, {output})",
    input=[
        q.pt_2,
        q.eta_2,
        q.phi_2,
        q.mass_2,
    ],
    output=[q.p4_2],
    scopes=["mm"],
)
LV_MM_reconstruction = Producer(
    name="LV_MM_reconstruction",
    call="lorentzvectors::CombineP4s({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2],
    output=[q.p4_mm],
    scopes=["mm"],
)
