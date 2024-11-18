from code_generation.producer import Producer
from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD

####################
# Set of producers used to Define Invariant mass of di-Muon event
####################

NumMuonCut = Producer(
    name="NumMuonCut",
    call="solution::CutNMuon({df}, {input}, {n_muon_req})",
    input=[nanoAOD.nMuon],
    output=[],
    scopes=["mm"],
)

CMuonCut = Producer(
    name="CMuonCut",
    call="solution::CutCMuons({df}, {input}, {C_muon_req})",
    input=[nanoAOD.Muon_charge],
    output=[],
    scopes=["mm"],
)

MuonInvMass = Producer(
    name="MuonInvMass",
    call="solution::GetInvariantMass({df}, {output}, {input})",
    input=[nanoAOD.Muon_pt, nanoAOD.Muon_eta, nanoAOD.Muon_phi, nanoAOD.Muon_mass],
    output=[q.Muon_InvMass],
    scopes=["mm"],
)

MuonCSum = Producer(
    name="MuonCSum",
    call="solution::MuonCSum({df}, {output}, {input})",
    input=[nanoAOD.Muon_charge],
    output=[q.MuonCSum],
    scopes=["mm"],
)