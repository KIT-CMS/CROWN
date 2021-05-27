import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup, VectorProducer

####################
# Set of producers used for selection of good taus
####################

TauPtCut = Producer(
    name="TauPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_tau_pt})",
    input=[nanoAOD.Tau_pt],
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
TauIDFilters = VectorProducer(
    name="TauIDFilters",
    call='physicsobject::tau::FilterTauID({df}, {output}, "{tau_id}", {tau_id_idx})',
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
    subproducers=[TauPtCut, TauEtaCut, TauDzCut, TauIDFilters],
)
