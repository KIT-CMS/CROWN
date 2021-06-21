import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import (
    Producer,
    ProducerGroup,
    VectorProducer,
    TriggerVectorProducer,
)

####################
# Set of producers used for trigger flags
####################

GenerateSingleMuonTriggerFlags = TriggerVectorProducer(
    name="GenerateSingleMuonTriggerFlags",
    call='trigger::GenerateSingleTriggerFlag({df}, {output}, {input}, "{hlt_path}", {ptcut}, {etacut}, {trigger_particle_id}, {filterbit}, {max_deltaR_triggermatch} )',
    input=[
        q.p4_1,
        nanoAOD.TriggerObject_bit,
        nanoAOD.TriggerObject_id,
        nanoAOD.TriggerObject_pt,
        nanoAOD.TriggerObject_eta,
        nanoAOD.TriggerObject_phi,
    ],
    output="flagname",
    scope=["mt"],
    vec_config="singlemoun_trigger",
)
GenerateCrossTriggerFlags = TriggerVectorProducer(
    name="GenerateCrossTriggerFlags",
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
    output="flagname",
    scope=["mt"],
    vec_config="cross_trigger",
)
MuTauTriggerFlags = ProducerGroup(
    name="MuTauTriggerFlags",
    call=None,
    input=None,
    output=None,
    scopes=["mt"],
    subproducers=[
        GenerateSingleMuonTriggerFlags,
        GenerateCrossTriggerFlags
    ],
)
