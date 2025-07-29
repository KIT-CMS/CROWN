from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for loosest selection of muons
####################

MuonPtCut = Producer(
    name="MuonPtCut",
    call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_muon_pt})",
    input=[nanoAOD.Muon_pt],
    output=[],
    scopes=["global"],
)
MuonEtaCut = Producer(
    name="MuonEtaCut",
    call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_muon_eta})",
    input=[nanoAOD.Muon_eta],
    output=[],
    scopes=["global"],
)
MuonDxyCut = Producer(
    name="MuonDxyCut",
    call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_muon_dxy})",
    input=[nanoAOD.Muon_dxy],
    output=[],
    scopes=["global"],
)
MuonDzCut = Producer(
    name="MuonDzCut",
    call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_muon_dz})",
    input=[nanoAOD.Muon_dz],
    output=[],
    scopes=["global"],
)
MuonIDCut = Producer(
    name="MuonIDCut",
    call="physicsobject::CutEqual<bool>({df}, {output}, {input}, true)",
    input=[nanoAOD.Muon_mediumId],
    output=[],
    scopes=["global"],
)
MuonIsoCut = Producer(
    name="MuonIsoCut",
    call="physicsobject::CutMax<float>({df}, {output}, {input}, {muon_iso_cut})",
    input=[nanoAOD.Muon_pfRelIso04_all],
    output=[],
    scopes=["global"],
)
BaseMuons = ProducerGroup(
    name="BaseMuons",
    call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
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
    call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_muon_pt})",
    input=[nanoAOD.Muon_pt],
    output=[],
    scopes=["em", "mt", "mm"],
)
GoodMuonEtaCut = Producer(
    name="GoodMuonEtaCut",
    call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_muon_eta})",
    input=[nanoAOD.Muon_eta],
    output=[],
    scopes=["em", "mt", "mm"],
)
GoodMuonIsoCut = Producer(
    name="GoodMuonIsoCut",
    call="physicsobject::CutMax<float>({df}, {output}, {input}, {muon_iso_cut})",
    input=[nanoAOD.Muon_pfRelIso04_all],
    output=[],
    scopes=["em", "mt", "mm"],
)
GoodMuons = ProducerGroup(
    name="GoodMuons",
    call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
    input=[q.base_muons_mask],
    output=[q.good_muons_mask],
    scopes=["em", "mt", "mm"],
    subproducers=[
        GoodMuonPtCut,
        GoodMuonEtaCut,
        GoodMuonIsoCut,
    ],
)
NumberOfGoodMuons = Producer(
    name="NumberOfGoodMuons",
    call="physicsobject::Count({df}, {output}, {input})",
    input=[q.good_muons_mask],
    output=[q.nmuons],
    scopes=["mt", "em", "mm"],
)
VetoMuons = Producer(
    name="VetoMuons",
    call="physicsobject::VetoSingleObject({df}, {output}, {input}, {muon_index_in_pair})",
    input=[q.base_muons_mask, q.dileptonpair],
    output=[q.veto_muons_mask],
    scopes=["em", "mt", "mm"],
)
VetoSecondMuon = Producer(
    name="VetoSecondMuon",
    call="physicsobject::VetoSingleObject({df}, {output}, {input}, {second_muon_index_in_pair})",
    input=[q.veto_muons_mask, q.dileptonpair],
    output=[q.veto_muons_mask_2],
    scopes=["mm"],
)

ExtraMuonsVeto = Producer(
    name="ExtraMuonsVeto",
    call="physicsobject::Veto({df}, {output}, {input})",
    input={
        "mm": [q.veto_muons_mask_2],
        "em": [q.veto_muons_mask],
        "et": [q.base_muons_mask],
        "mt": [q.veto_muons_mask],
        "tt": [q.base_muons_mask],
    },
    output=[q.muon_veto_flag],
    scopes=["em", "et", "mt", "tt", "mm"],
)

####################
# Set of producers used for di-muon veto
####################

DiMuonVetoPtCut = Producer(
    name="DiMuonVetoPtCut",
    call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_dimuonveto_pt})",
    input=[nanoAOD.Muon_pt],
    output=[],
    scopes=["global"],
)
DiMuonVetoIDCut = Producer(
    name="DiMuonVetoIDCut",
    call="physicsobject::CutEqual<bool>({df}, {output}, {input}, true)",
    input=[nanoAOD.Muon_looseId],
    output=[],
    scopes=["global"],
)
DiMuonVetoMuons = ProducerGroup(
    name="DiMuonVetoMuons",
    call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
    input=MuonEtaCut.output + MuonDxyCut.output + MuonDzCut.output + MuonIsoCut.output,
    output=[],
    scopes=["global"],
    subproducers=[
        DiMuonVetoPtCut,
        DiMuonVetoIDCut,
    ],
)
DiMuonVeto = ProducerGroup(
    name="DiMuonVeto",
    call="physicsobject::LeptonPairVeto({df}, {output}, {input}, {dileptonveto_dR})",
    input=[
        nanoAOD.Muon_pt,
        nanoAOD.Muon_eta,
        nanoAOD.Muon_phi,
        nanoAOD.Muon_mass,
        nanoAOD.Muon_charge,
    ],
    output=[q.dimuon_veto],
    scopes=["global"],
    subproducers=[DiMuonVetoMuons],
)
