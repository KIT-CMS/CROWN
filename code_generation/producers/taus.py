import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup


#############
# Tau ID cuts
#############

VsJetTauIDCut = Producer(
    name="VsJetTauIDCut",
    call="physicsobject::tau::CutTauID({df}, {output}, {input}, {vsjet_tau_id_bit})",
    input=[nanoAOD.Tau_ID_vsJet],
    output=[],
    scopes=["et", "mt", "tt"],
)
VsElectronTauIDCut = Producer(
    name="VsElectronTauIDCut",
    call="physicsobject::tau::CutTauID({df}, {output}, {input}, {vsele_tau_id_bit})",
    input=[nanoAOD.Tau_ID_vsEle],
    output=[],
    scopes=["et", "mt", "tt"],
)
VsMuonTauIDCut = Producer(
    name="VsMuonTauIDCut",
    call="physicsobject::tau::CutTauID({df}, {output}, {input}, {vsmu_tau_id_bit})",
    input=[nanoAOD.Tau_ID_vsMu],
    output=[],
    scopes=["et", "mt", "tt"],
)

####################
# Set of producers used for selection of good taus
####################

TauPtCorrection_byValue = Producer(
    name="TauPtCorrection",
    call="physicsobject::tau::PtCorrection_byValue({df}, {output}, {input}, {tau_ES_shift_DM0}, {tau_ES_shift_DM1}, {tau_ES_shift_DM10}, {tau_ES_shift_DM11})",
    input=[
        nanoAOD.Tau_pt,
        nanoAOD.Tau_decayMode,
    ],
    output=[q.Tau_pt_corrected],
    scopes=["global"],
)
TauPtCorrection = Producer(
    name="TauPtCorrection",
    call='physicsobject::tau::PtCorrection({df}, {output}, {input}, "{tau_sf_file}", "{tau_ES_json_name}", "{tau_id_algorithm}", "{tau_ES_shift_DM0}", "{tau_ES_shift_DM1}", "{tau_ES_shift_DM10}", "{tau_ES_shift_DM11}", {vec_open}{tau_dms}{vec_close})',
    input=[
        nanoAOD.Tau_pt,
        nanoAOD.Tau_eta,
        nanoAOD.Tau_decayMode,
        nanoAOD.Tau_genMatch,
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
TauEnergyCorrection_byValue = ProducerGroup(
    name="TauEnergyCorrection",
    call=None,
    input=None,
    output=None,
    scopes=["global"],
    subproducers=[TauPtCorrection_byValue, TauMassCorrection],
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
# TauIDCuts = VectorProducer(
#     name="TauIDCuts",
#     call='physicsobject::tau::CutTauID({df}, {output}, "{tau_id}", {tau_id_idx})',
#     input=[],
#     output=[],
#     scopes=["global"],
#     vec_configs=["tau_id", "tau_id_idx"],
# )
BaseTaus = ProducerGroup(
    name="BaseTaus",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.base_taus_mask],
    scopes=["global"],
    subproducers=[TauPtCut, TauEtaCut, TauDzCut, TauDMCut],
)
######
# Good taus selection for nonglobal scope
######
GoodTauPtCut = Producer(
    name="GoodTauPtCut",
    call="physicsobject::CutPt({df}, {input}, {output}, {min_tau_pt})",
    input=[q.Tau_pt_corrected],
    output=[],
    scopes=["et", "mt", "tt"],
)
GoodTauEtaCut = Producer(
    name="GoodTauEtaCut",
    call="physicsobject::CutEta({df}, {input}, {output}, {max_tau_eta})",
    input=[nanoAOD.Tau_eta],
    output=[],
    scopes=["et", "mt", "tt"],
)
GoodTauDzCut = Producer(
    name="GoodTauDzCut",
    call="physicsobject::CutDz({df}, {input}, {output}, {max_tau_dz})",
    input=[nanoAOD.Tau_dz],
    output=[],
    scopes=["et", "mt", "tt"],
)
GoodTauDMCut = Producer(
    name="GoodTauDMCut",
    call="physicsobject::tau::CutDecayModes({df}, {output}, {input}, {vec_open}{tau_dms}{vec_close})",
    input=[nanoAOD.Tau_decayMode],
    output=[],
    scopes=["et", "mt", "tt"],
)
GoodTaus = ProducerGroup(
    name="GoodTaus",
    call="physicsobject::CombineMasks({df}, {output}, {input})",
    input=[],
    output=[q.good_taus_mask],
    scopes=["et", "mt", "tt"],
    subproducers=[
        GoodTauPtCut,
        GoodTauEtaCut,
        GoodTauDzCut,
        GoodTauDMCut,
        VsJetTauIDCut,
        VsElectronTauIDCut,
        VsMuonTauIDCut,
    ],
)
NumberOfGoodTaus = Producer(
    name="NumberOfGoodTaus",
    call="quantities::NumberOfGoodLeptons({df}, {output}, {input})",
    input=[q.good_taus_mask],
    output=[q.ntaus],
    scopes=["mt", "et", "tt"],
)
