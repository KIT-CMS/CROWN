from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of general producers for DiTauPair Quantities
####################

pt_1 = Producer(
    name="pt_1",
    call="lorentzvector::GetPt({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.pt_1],
    scopes=["mm"],
)
pt_2 = Producer(
    name="pt_2",
    call="lorentzvector::GetPt({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.pt_2],
    scopes=["mm"],
)
eta_1 = Producer(
    name="eta_1",
    call="lorentzvector::GetEta({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.eta_1],
    scopes=["mm"],
)
eta_2 = Producer(
    name="eta_2",
    call="lorentzvector::GetEta({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.eta_2],
    scopes=["mm"],
)
phi_1 = Producer(
    name="phi_1",
    call="lorentzvector::GetPhi({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.phi_1],
    scopes=["mm"],
)
phi_2 = Producer(
    name="phi_2",
    call="lorentzvector::GetPhi({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.phi_2],
    scopes=["mm"],
)
mass_1 = Producer(
    name="mass_1",
    call="lorentzvector::GetMass({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.mass_1],
    scopes=["mm"],
)
mass_2 = Producer(
    name="mass_2",
    call="lorentzvector::GetMass({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.mass_2],
    scopes=["mm"],
)
muon_q_1 = Producer(
    name="muon_q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_charge],
    output=[q.q_1],
    scopes=["mm"],
)
muon_q_2 = Producer(
    name="muon_q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_charge],
    output=[q.q_2],
    scopes=["mm"],
)

UnrollMuLV1 = ProducerGroup(
    name="UnrollMuLV1",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[
        pt_1,
        eta_1,
        phi_1,
        mass_1,
        muon_q_1,
    ],
)
UnrollMuLV2 = ProducerGroup(
    name="UnrollMuLV2",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[
        pt_2,
        eta_2,
        phi_2,
        mass_2,
        muon_q_2,
    ],
)

mm_pair_mass = Producer(
    name="mm_pair_mass",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.mm_pair_mass],
    scopes=["mm"],
)
mm_pair_pt = Producer(
    name="mm_pair_pt",
    call="quantities::pt_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.mm_pair_pt],
    scopes=["mm"],
)

#####################
# Producer Groups
#####################

MMDiTauPairQuantities = ProducerGroup(
    name="MMDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[UnrollMuLV1, UnrollMuLV2, mm_pair_mass, mm_pair_pt],
)
