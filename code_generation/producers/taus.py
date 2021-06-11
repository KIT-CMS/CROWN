import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup, VectorProducer

####################
# Set of producers used for selection of good taus
####################

TauPtCorrection = Producer(
    name="TauPtCorrection",
    call="physicsobject::tau::PtCorrection({df}, {output}, {input}, {tau_ES_shift_DM0}, {tau_ES_shift_DM1}, {tau_ES_shift_DM10}, {tau_ES_shift_DM11})",
    input=[
        nanoAOD.Tau_pt,
        nanoAOD.Tau_decayMode,
    ],
    output=[q.Tau_pt_corrected],
    scopes=["global"],
)
TauMassCorrection = Producer(
    name="TauMassCorrection",
    call="physicsobject::ObjectMassCorrectionWithPt({df}, {output}, {input})",
    input=[
        nanoAOD.Tau_mass,
        nanoAOD.Tau_pt,
        q.Tau_pt_corrected,
    ],
    output=[q.Tau_mass_corrected],
    scopes=["global"],
)
TauEnergyCorrection = ProducerGroup(
    name="TauEnergyCorrection",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[TauPtCorrection, TauMassCorrection],
)
TauPtCut = Producer(
    name="TauPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_tau_pt})",
    input=[q.Tau_pt_corrected],
    output=[],
    scopes=["global"],
)
TauEtaCut = Producer(
    name="TauEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_tau_eta})",
    input=[nanoAOD.Tau_eta],
    output=[],
    scopes=["global"],
)
TauDzCut = Producer(
    name="TauDzCut",
    call="physicsobject::CutDz({df}, {input}, {output}, {max_tau_dz})",
    input=[nanoAOD.Tau_dz],
    output=[],
    scopes=["global"],
)
TauDMCut = Producer(
    name="TauDMCut",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{tau_dms}{vec_close})",
    input=[nanoAOD.Tau_decayMode],
    output=[],
    scopes=["global"],
)
TauIDCuts = VectorProducer(
    name="TauIDCuts",
    call='physicsobject::tau::CutTauID({df}, {output}, "{tau_id}", {tau_id_idx})',
    input=[],
    output=[],
    scopes=["global"],
    vec_configs=["tau_id", "tau_id_idx"],
)
GoodTaus = ProducerGroup(
    name="GoodTaus",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_taus_mask],
    scopes=["global"],
    subproducers=[TauPtCut, TauEtaCut, TauDzCut, TauDMCut, TauIDCuts],
)
