from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for loosest selection of muons
####################

MuonPtCut = Producer(
    name="MuonPtCut",
    call="physicsobject::CutMin<float>({df}, {output}, {input}, {muon_min_pt})",
    input=[nanoAOD.Muon_pt],
    output=[],
    scopes=["global"],
)
MuonEtaCut = Producer(
    name="MuonEtaCut",
    call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {muon_max_eta})",
    input=[nanoAOD.Muon_eta],
    output=[],
    scopes=["global"],
)
MuonDxyCut = Producer(
    name="MuonDxyCut",
    call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {muon_max_dxy})",
    input=[nanoAOD.Muon_dxy],
    output=[],
    scopes=["global"],
)
MuonDzCut = Producer(
    name="MuonDzCut",
    call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {muon_max_dz})",
    input=[nanoAOD.Muon_dz],
    output=[],
    scopes=["global"],
)
MuonIDCut = Producer(
    name="MuonIDCut",
    call='physicsobject::CutEqual<bool>({df}, {output}, {input}, true)',
    input=[nanoAOD.Muon_mediumId],
    output=[],
    scopes=["global"],
)
MuonIsoCut = Producer(
    name="MuonIsoCut",
    call="physicsobject::CutMax<float>({df}, {output}, {input}, {muon_max_iso})",
    input=[nanoAOD.Muon_iso],
    output=[],
    scopes=["global"],
)
BaseMuons = ProducerGroup(
    name="BaseMuons",
    call='physicsobject::CombineMasks({df}, {output}, {input}, "all")',
    input=[],
    output=[q.base_muons_mask],
    scopes=["global"],
    subproducers=[
        MuonPtCut,
        MuonEtaCut,
        MuonDxyCut,
        MuonDzCut,
        MuonIDCut,
        MuonIsoCut,
    ],
)

####################
# Set of producers used for more specific selection of muons in channels
####################

GoodMuonPtCut = Producer(
    name="GoodMuonPtCut",
    call="physicsobject::CutMin<float>({df}, {output}, {input}, {muon_min_pt})",
    input=[nanoAOD.Muon_pt],
    output=[],
    scopes=["mm"],
)
GoodMuonEtaCut = Producer(
    name="GoodMuonEtaCut",
    call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {muon_max_eta})",
    input=[nanoAOD.Muon_eta],
    output=[],
    scopes=["mm"],
)
GoodMuonIsoCut = Producer(
    name="GoodMuonIsoCut",
    call="physicsobject::CutMax<float>({df}, {output}, {input}, {muon_max_iso})",
    input=[nanoAOD.Muon_iso],
    output=[],
    scopes=["mm"],
)
GoodMuons = ProducerGroup(
    name="GoodMuons",
    call='physicsobject::CombineMasks({df}, {output}, {input}, "all")',
    input=[q.base_muons_mask],
    output=[q.good_muons_mask],
    scopes=["mm"],
    subproducers=[
        GoodMuonPtCut,
        GoodMuonEtaCut,
        GoodMuonIsoCut,
    ],
)

NumberOfGoodMuons = Producer(
    name="NumberOfGoodMuons",
    call="quantities::NumberOfGoodLeptons({df}, {output}, {input})",
    input=[q.good_muons_mask],
    output=[q.nmuons],
    scopes=["mm"],
)
