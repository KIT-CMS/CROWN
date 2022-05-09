import code_generation.quantities.output as q
import code_generation.quantities.tagandprobe_output as tp_q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producers.muons import (
    MuonPtCut,
    MuonEtaCut,
    MuonDxyCut,
    MuonDzCut,
    MuonIsoCut,
)
from code_generation.producer import (
    Producer,
    Filter,
    ProducerGroup,
    TriggerVectorProducer,
)

BaseMuons = ProducerGroup(
    name="BaseMuons",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.base_muons_mask],
    scopes=["global"],
    subproducers=[
        MuonPtCut,
        MuonEtaCut,
        MuonDxyCut,
        MuonDzCut,
        MuonIsoCut,
    ],
)

MMTagAndProbePairSelection = Producer(
    name="MMPairSelection",
    call="pairselection::mumu::TagAndProbePairSelection({df}, {input_vec}, {output})",
    input=[
        nanoAOD.Muon_pt,
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

MMTagAndProbePairs = ProducerGroup(
    name="UnrollMuLV1",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[
        MMTagAndProbePairSelection,
        GoodMMPairFlag,
        GoodMMPairFilter,
    ],
)
MMSingleMuonTriggerFlags_1 = TriggerVectorProducer(
    name="MMGenerateSingleMuonTriggerFlags",
    call='trigger::GenerateSingleTriggerFlag({df}, {output}, {input}, "{hlt_path}", {ptcut}, {etacut}, {trigger_particle_id}, {filterbit}, {max_deltaR_triggermatch} )',
    input=[
        q.p4_1,
        nanoAOD.TriggerObject_bit,
        nanoAOD.TriggerObject_id,
        nanoAOD.TriggerObject_pt,
        nanoAOD.TriggerObject_eta,
        nanoAOD.TriggerObject_phi,
    ],
    output="flagname_1",
    scope=["mm"],
    vec_config="singlemoun_trigger",
)
MMSingleMuonTriggerFlags_2 = TriggerVectorProducer(
    name="MMGenerateSingleMuonTriggerFlags",
    call='trigger::GenerateSingleTriggerFlag({df}, {output}, {input}, "{hlt_path}", {ptcut}, {etacut}, {trigger_particle_id}, {filterbit}, {max_deltaR_triggermatch} )',
    input=[
        q.p4_2,
        nanoAOD.TriggerObject_bit,
        nanoAOD.TriggerObject_id,
        nanoAOD.TriggerObject_pt,
        nanoAOD.TriggerObject_eta,
        nanoAOD.TriggerObject_phi,
    ],
    output="flagname_2",
    scope=["mm"],
    vec_config="singlemoun_trigger",
)
MMDoubleMuonTriggerFlags_1 = TriggerVectorProducer(
    name="MMDoubleMuonTriggerFlags_1",
    call='trigger::GenerateDoubleTriggerFlag({df}, {output}, {input}, "{hlt_path}", {p1_ptcut}, {p2_ptcut}, {p1_etacut}, {p2_etacut}, {p1_trigger_particle_id}, {p2_trigger_particle_id}, {p1_filterbit}, {p2_filterbit}, {max_deltaR_triggermatch})',
    input=[
        q.p4_1,
        q.p4_2,
        nanoAOD.TriggerObject_bit,
        nanoAOD.TriggerObject_id,
        nanoAOD.TriggerObject_pt,
        nanoAOD.TriggerObject_eta,
        nanoAOD.TriggerObject_phi,
    ],
    output="flagname_1",
    scope=["mm"],
    vec_config="doublemuon_trigger",
)
MMDoubleMuonTriggerFlags_2 = TriggerVectorProducer(
    name="MMDoubleMuonTriggerFlags_2",
    call='trigger::GenerateDoubleTriggerFlag({df}, {output}, {input}, "{hlt_path}", {p1_ptcut}, {p2_ptcut}, {p1_etacut}, {p2_etacut}, {p1_trigger_particle_id}, {p2_trigger_particle_id}, {p1_filterbit}, {p2_filterbit}, {max_deltaR_triggermatch})',
    input=[
        q.p4_2,
        q.p4_1,
        nanoAOD.TriggerObject_bit,
        nanoAOD.TriggerObject_id,
        nanoAOD.TriggerObject_pt,
        nanoAOD.TriggerObject_eta,
        nanoAOD.TriggerObject_phi,
    ],
    output="flagname_2",
    scope=["mm"],
    vec_config="doublemuon_trigger",
)

## Producers to writeout the id variables for the tag and probe pairs
MuonID_Loose_1 = Producer(
    name="MuonID_Loose_1",
    call="quantities::muon::id({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_id_loose],
    output=[tp_q.id_loose_1],
    scopes=["mm"],
)
MuonID_Loose_2 = Producer(
    name="MuonID_Loose_2",
    call="quantities::muon::id({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Muon_id_loose],
    output=[tp_q.id_loose_2],
    scopes=["mm"],
)
MuonID_Medium_1 = Producer(
    name="MuonID_Medium_1",
    call="quantities::muon::id({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_id_medium],
    output=[tp_q.id_medium_1],
    scopes=["mm"],
)
MuonID_Medium_2 = Producer(
    name="MuonID_Medium_2",
    call="quantities::muon::id({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Muon_id_medium],
    output=[tp_q.id_medium_2],
    scopes=["mm"],
)
MuonID_Tight_1 = Producer(
    name="MuonID_Tight_1",
    call="quantities::muon::id({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_id_tight],
    output=[tp_q.id_tight_1],
    scopes=["mm"],
)
MuonID_Tight_2 = Producer(
    name="MuonID_Tight_2",
    call="quantities::muon::id({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Muon_id_tight],
    output=[tp_q.id_tight_2],
    scopes=["mm"],
)
MuonIDs = ProducerGroup(
    name="MuonIDs",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[
        MuonID_Loose_1,
        MuonID_Loose_2,
        MuonID_Medium_1,
        MuonID_Medium_2,
        MuonID_Tight_1,
        MuonID_Tight_2,
    ],
)
