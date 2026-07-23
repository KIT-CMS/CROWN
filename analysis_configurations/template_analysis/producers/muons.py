from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

####################
# Set of producers used for loosest selection of muons
####################

with defaults(scopes=["global"]):
    with defaults(output=[]):
        MuonPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {muon_min_pt})",
            input=[nanoAOD.Muon_pt],
        )
        MuonEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {muon_max_eta})",
            input=[nanoAOD.Muon_eta],
        )
        MuonDxyCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {muon_max_dxy})",
            input=[nanoAOD.Muon_dxy],
        )
        MuonDzCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {muon_max_dz})",
            input=[nanoAOD.Muon_dz],
        )
        MuonIDCut = Producer(
            call="physicsobject::CutEqual<bool>({df}, {output}, {input}, true)",
            input=[nanoAOD.Muon_mediumId],
        )
        MuonIsoCut = Producer(
            call="physicsobject::CutMax<float>({df}, {output}, {input}, {muon_max_iso})",
            input=[nanoAOD.Muon_pfRelIso04_all],
        )
    BaseMuons = ProducerGroup(
        call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
        input=[],
        output=[q.base_muons_mask],
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

with defaults(scopes=["mm"]):
    with defaults(output=[]):
        GoodMuonPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {muon_min_pt})",
            input=[nanoAOD.Muon_pt],
        )
        GoodMuonEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {muon_max_eta})",
            input=[nanoAOD.Muon_eta],
        )
        GoodMuonIsoCut = Producer(
            call="physicsobject::CutMax<float>({df}, {output}, {input}, {muon_max_iso})",
            input=[nanoAOD.Muon_pfRelIso04_all],
        )
    GoodMuons = ProducerGroup(
        call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
        input=[q.base_muons_mask],
        output=[q.good_muons_mask],
        subproducers=[
            GoodMuonPtCut,
            GoodMuonEtaCut,
            GoodMuonIsoCut,
        ],
    )
    NumberOfGoodMuons = Producer(
        call="physicsobject::Count({df}, {output}, {input})",
        input=[q.good_muons_mask],
        output=[q.nmuons],
    )
