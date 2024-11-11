from code_generation.producer import Producer, ProducerGroup
from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD

####################
# Set of producers used to Define Invariant mass of di-Muon event
####################

NumMuonCut = Producer(
    name="NumMuonCut",
    call="example::CutNMuon({df}, {input}, {n_muon_req})",
    input=[nanoAOD.nMuon],
    output=[],
    scopes=["leptonic"],
)

CMuonCut = Producer(
    name="CMuonCut",
    call="example::CutCMuons({df}, {input}, {C_muon_req})",
    input=[nanoAOD.Muon_charge],
    output=[],
    scopes=["leptonic"],
)

MuonInvMass = Producer(
    name="MuonInvMass",
    call="example::GetInvariantMass({df}, {output}, {input})",
    input=[nanoAOD.Muon_pt, nanoAOD.Muon_eta, nanoAOD.Muon_phi, nanoAOD.Muon_mass],
    output=[q.Muon_InvMass],
    scopes=["leptonic"],
)