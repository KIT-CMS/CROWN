from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

with defaults(scopes=["global"], output=[]):
    TauPtCut = Producer(
        call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_tau_pt})",
        input=[q.Tau_pt_corrected],
    )
    TauEtaCut = Producer(
        call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_tau_eta})",
        input=[nanoAOD.Tau_eta],
    )
    TauDzCut = Producer(
        call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_tau_dz})",
        input=[nanoAOD.Tau_dz],
    )
    TauDMCut = Producer(
        call="physicsobject::CutQuantity<int>({df}, {output}, {input}, {vec_open}{tau_dms}{vec_close})",
        input=[nanoAOD.Tau_decayMode],
    )

#############
# Tau ID cuts
#############
with defaults(scopes=["et", "mt", "tt"]):
    with defaults(output=[]):
        VsJetTauIDCut = Producer(
            call="physicsobject::CutBitmask({df}, {output}, {input}, {vsjet_tau_id_bit})",
            input=[nanoAOD.Tau_idDeepTau2017v2p1VSjet],
        )
        VsElectronTauIDCut = Producer(
            call="physicsobject::CutBitmask({df}, {output}, {input}, {vsele_tau_id_bit})",
            input=[nanoAOD.Tau_idDeepTau2017v2p1VSe],
        )
        VsMuonTauIDCut = Producer(
            call="physicsobject::CutBitmask({df}, {output}, {input}, {vsmu_tau_id_bit})",
            input=[nanoAOD.Tau_idDeepTau2017v2p1VSmu],
        )

    ####################
    # Set of producers used for selection of good taus
    ####################
    TauPtCorrection_byValue = Producer(
        call="""embedding::tau::PtCorrection_byValue(
            {df}, 
            {output}, 
            {input}, 
            {tau_ES_shift_DM0}, 
            {tau_ES_shift_DM1}, 
            {tau_ES_shift_DM10}, 
            {tau_ES_shift_DM11})
            """,
        input=[
            nanoAOD.Tau_pt,
            nanoAOD.Tau_decayMode,
        ],
        output=[q.Tau_pt_corrected],
    )
    TauPtCorrection_eleFake = Producer(
        call="""physicsobject::tau::PtCorrectionMC_eleFake(
            {df}, 
            correctionManager, 
            {output}, 
            {input}, 
            "{tau_sf_file}", 
            "{tau_ES_json_name}", 
            "{tau_id_algorithm}", 
            "{tau_elefake_es_DM0_barrel}", 
            "{tau_elefake_es_DM1_barrel}", 
            "{tau_elefake_es_DM0_endcap}", 
            "{tau_elefake_es_DM1_endcap}")
            """,
        input=[
            nanoAOD.Tau_pt,
            nanoAOD.Tau_eta,
            nanoAOD.Tau_decayMode,
            nanoAOD.Tau_genPartFlav,
        ],
        output=[q.Tau_pt_ele_corrected],
    )
    TauPtCorrection_muFake = Producer(
        call="""physicsobject::tau::PtCorrectionMC_muFake(
            {df}, 
            correctionManager, 
            {output}, 
            {input}, 
            "{tau_sf_file}", 
            "{tau_ES_json_name}", 
            "{tau_id_algorithm}", 
            "{tau_mufake_es}")
            """,
        input=[
            q.Tau_pt_ele_corrected,
            nanoAOD.Tau_eta,
            nanoAOD.Tau_decayMode,
            nanoAOD.Tau_genPartFlav,
        ],
        output=[q.Tau_pt_ele_mu_corrected],
    )
    TauPtCorrection_genTau = Producer(
        call="""physicsobject::tau::PtCorrectionMC_genuineTau(
            {df}, 
            correctionManager, 
            {output}, 
            {input}, 
            "{tau_sf_file}", 
            "{tau_ES_json_name}", 
            "{tau_id_algorithm}", 
            "{tau_ES_shift_DM0}", 
            "{tau_ES_shift_DM1}", 
            "{tau_ES_shift_DM10}", 
            "{tau_ES_shift_DM11}")
            """,
        input=[
            q.Tau_pt_ele_mu_corrected,
            nanoAOD.Tau_eta,
            nanoAOD.Tau_decayMode,
            nanoAOD.Tau_genPartFlav,
        ],
        output=[q.Tau_pt_corrected],
    )
    with defaults(
        call="event::quantity::Rename<ROOT::RVec<float>>({df}, {output}, {input})"
    ):
        TauPtCorrection_data = Producer(
            input=[nanoAOD.Tau_pt],
            output=[q.Tau_pt_corrected],
        )
        TauMassCorrection_data = Producer(
            input=[nanoAOD.Tau_mass],
            output=[q.Tau_mass_corrected],
        )
    TauMassCorrection = Producer(
        call="physicsobject::MassCorrectionWithPt({df}, {output}, {input})",
        input=[
            nanoAOD.Tau_mass,
            nanoAOD.Tau_pt,
            q.Tau_pt_corrected,
        ],
        output=[q.Tau_mass_corrected],
    )
    with defaults(call=None, input=None, output=None):
        TauEnergyCorrection_byValue = ProducerGroup(
            subproducers=[
                TauPtCorrection_eleFake,
                TauPtCorrection_byValue,
                TauMassCorrection,
            ],
        )
        TauEnergyCorrection = ProducerGroup(
            subproducers=[
                TauPtCorrection_eleFake,
                TauPtCorrection_muFake,
                TauPtCorrection_genTau,
                TauMassCorrection,
            ],
        )
        TauEnergyCorrection_data = ProducerGroup(
            subproducers=[TauPtCorrection_data, TauMassCorrection_data],
        )

    ######
    # Good taus selection for nonglobal scope
    ######
    with defaults(output=[]):
        GoodTauPtCut = Producer(
            call="physicsobject::CutMin<float>({df}, {output}, {input}, {min_tau_pt})",
            input=[q.Tau_pt_corrected],
        )
        GoodTauEtaCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_tau_eta})",
            input=[nanoAOD.Tau_eta],
        )
        GoodTauDzCut = Producer(
            call="physicsobject::CutAbsMax<float>({df}, {output}, {input}, {max_tau_dz})",
            input=[nanoAOD.Tau_dz],
        )
        GoodTauDMCut = Producer(
            call="physicsobject::CutQuantity<int>({df}, {output}, {input}, {vec_open}{tau_dms}{vec_close})",
            input=[nanoAOD.Tau_decayMode],
        )
    GoodTaus = ProducerGroup(
        call='physicsobject::CombineMasks({df}, {output}, {input}, "all_of")',
        input=[],
        output=[q.good_taus_mask],
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
        call="physicsobject::Count({df}, {output}, {input})",
        input=[q.good_taus_mask],
        output=[q.ntaus],
    )
