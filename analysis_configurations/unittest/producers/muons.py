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
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_muon_pt})",
            input=[nanoAOD.Muon_pt],
        )
        MuonEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_muon_eta})",
            input=[nanoAOD.Muon_eta],
        )
        MuonDxyCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_muon_dxy})",
            input=[nanoAOD.Muon_dxy],
        )
        MuonDzCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_muon_dz})",
            input=[nanoAOD.Muon_dz],
        )
        MuonIDCut = Producer(
            call="physicsobject::CutEqual<bool>({df}, {output}, {input}, true)",
            input=[nanoAOD.Muon_mediumId],
        )
        MuonIsoCut = Producer(
            call="physicsobject::CutMax<float>({df}, {output}, {input}, {muon_iso_cut})",
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

with defaults(scopes=["em", "mt", "mm"]):
    with defaults(output=[]):
        GoodMuonPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_muon_pt})",
            input=[nanoAOD.Muon_pt],
        )
        GoodMuonEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_muon_eta})",
            input=[nanoAOD.Muon_eta],
        )
        GoodMuonIsoCut = Producer(
            call="physicsobject::CutMax<float>({df}, {output}, {input}, {muon_iso_cut})",
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
    VetoMuons = Producer(
        call="physicsobject::VetoSingleObject({df}, {output}, {input}, {muon_index_in_pair})",
        input=[q.base_muons_mask, q.dileptonpair],
        output=[q.veto_muons_mask],
    )
VetoSecondMuon = Producer(
    call="physicsobject::VetoSingleObject({df}, {output}, {input}, {second_muon_index_in_pair})",
    input=[q.veto_muons_mask, q.dileptonpair],
    output=[q.veto_muons_mask_2],
    scopes=["mm"],
)
ExtraMuonsVeto = Producer(
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

with defaults(scopes=["global"]):
    with defaults(output=[]):
        DiMuonVetoPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_dimuonveto_pt})",
            input=[nanoAOD.Muon_pt],
        )
        DiMuonVetoIDCut = Producer(
            call="physicsobject::CutEqual<bool>({df}, {output}, {input}, true)",
            input=[nanoAOD.Muon_looseId],
        )
        DiMuonVetoMuons = ProducerGroup(
            call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
            input=MuonEtaCut.output
            + MuonDxyCut.output
            + MuonDzCut.output
            + MuonIsoCut.output,
            subproducers=[
                DiMuonVetoPtCut,
                DiMuonVetoIDCut,
            ],
        )
    DiMuonVeto = ProducerGroup(
        call="physicsobject::LeptonPairVeto({df}, {output}, {input}, {dileptonveto_dR})",
        input=[
            nanoAOD.Muon_pt,
            nanoAOD.Muon_eta,
            nanoAOD.Muon_phi,
            nanoAOD.Muon_mass,
            nanoAOD.Muon_charge,
        ],
        output=[q.dimuon_veto],
        subproducers=[DiMuonVetoMuons],
    )
