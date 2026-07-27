from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import ExtendedVectorProducer
from code_generation.helpers import defaults

####################
# Set of producers used for trigger flags
####################

_nanoAOD_TrigObj_vars = [
    nanoAOD.TrigObj_pt,
    nanoAOD.TrigObj_eta,
    nanoAOD.TrigObj_phi,
    nanoAOD.TrigObj_id,
    nanoAOD.TrigObj_filterBits,
]

with defaults(
    call="""trigger::SingleObjectFlag(
        {df},
        {output},
        {input},
        "{hlt_path}",
        {ptcut},
        {etacut},
        {trigger_particle_id},
        {vec_open}{filterbit}{vec_close},
        {max_deltaR_triggermatch})
    """,
    output="flagname",
):
    MMGenerateSingleMuonTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1] + _nanoAOD_TrigObj_vars,
        scopes=["mm"],
        vec_config="singlemoun_trigger",
    )
    MTGenerateSingleMuonTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1] + _nanoAOD_TrigObj_vars,
        scopes=["mt"],
        vec_config="singlemoun_trigger",
    )
    ETGenerateSingleElectronTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1] + _nanoAOD_TrigObj_vars,
        scopes=["et"],
        vec_config="singleelectron_trigger",
    )
    EMGenerateSingleElectronTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1] + _nanoAOD_TrigObj_vars,
        scopes=["em"],
        vec_config="singleelectron_trigger",
    )
    GenerateSingleLeadingTauTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1] + _nanoAOD_TrigObj_vars,
        scopes=["tt"],
        vec_config="singletau_trigger_leading",
    )
    GenerateSingleTrailingTauTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_2] + _nanoAOD_TrigObj_vars,
        scopes=["et", "mt", "tt"],
        vec_config="singletau_trigger_trailing",
    )
    EMGenerateSingleMuonTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_2] + _nanoAOD_TrigObj_vars,
        scopes=["em"],
        vec_config="singlemoun_trigger",
    )
with defaults(
    call="""trigger::DoubleObjectFlag(
        {df},
        {output},
        {input},
        "{hlt_path}",
        {p1_ptcut},
        {p2_ptcut},
        {p1_etacut},
        {p2_etacut},
        {p1_trigger_particle_id},
        {p2_trigger_particle_id},
        {vec_open}{p1_filterbit}{vec_close},
        {vec_open}{p2_filterbit}{vec_close},
        {max_deltaR_triggermatch})
    """,
    output="flagname",
):
    MTGenerateCrossTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1, q.p4_2] + _nanoAOD_TrigObj_vars,
        scopes=["mt"],
        vec_config="mutau_cross_trigger",
    )
    ETGenerateCrossTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1, q.p4_2] + _nanoAOD_TrigObj_vars,
        scopes=["et"],
        vec_config="eltau_cross_trigger",
    )
    TTGenerateDoubleTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1, q.p4_2] + _nanoAOD_TrigObj_vars,
        scopes=["tt"],
        vec_config="doubletau_trigger",
    )
    EMGenerateCrossTriggerFlags = ExtendedVectorProducer(
        input=[q.p4_1, q.p4_2] + _nanoAOD_TrigObj_vars,
        scopes=["em"],
        vec_config="elmu_cross_trigger",
    )
