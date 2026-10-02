from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

####################
# Set of producers used for contruction of met related quantities
####################
with defaults(scopes=["global"]):
    with defaults(call="lorentzvector::BuildMET({df}, {output}, {input})"):
        BuildMetVector = Producer(
            input=[
                nanoAOD.PuppiMET_pt,
                nanoAOD.PuppiMET_phi,
            ],
            output=[q.met_p4],
        )
        BuildPFMetVector = Producer(
            input=[
                nanoAOD.MET_pt,
                nanoAOD.MET_phi,
            ],
            output=[q.pfmet_p4],
        )
    with defaults(call="event::quantity::Rename<float>({df}, {output}, {input})"):
        MetCov00 = Producer(
            input=[nanoAOD.MET_covXX],
            output=[q.metcov00],
        )
        MetCov01 = Producer(
            input=[nanoAOD.MET_covXY],
            output=[q.metcov01],
        )
        MetCov10 = Producer(
            input=[nanoAOD.MET_covXY],
            output=[q.metcov10],
        )
        MetCov11 = Producer(
            input=[nanoAOD.MET_covYY],
            output=[q.metcov11],
        )
        MetSumEt = Producer(
            input=[nanoAOD.PuppiMET_sumEt],
            output=[q.metSumEt],
        )
    with defaults(input=[q.met_p4]):
        MetPt_uncorrected = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.met_uncorrected],
        )
        MetPhi_uncorrected = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.metphi_uncorrected],
        )
    with defaults(input=[q.pfmet_p4]):
        PFMetPt_uncorrected = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.pfmet_uncorrected],
        )
        PFMetPhi_uncorrected = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.pfmetphi_uncorrected],
        )
    with defaults(
        input=[
            nanoAOD.GenPart_pt,
            nanoAOD.GenPart_eta,
            nanoAOD.GenPart_phi,
            nanoAOD.GenPart_mass,
            nanoAOD.GenPart_pdgId,
            nanoAOD.GenPart_status,
            nanoAOD.GenPart_statusFlags,
        ]
    ):
        CalculateGenBosonVector = Producer(
            call="genparticles::GetBoson({df}, {output}, {input}, {is_data})",
            output=[q.recoil_genboson_p4],
        )
        CalculateVisGenBosonVector = Producer(
            call="genparticles::GetVisibleBoson({df}, {output}, {input}, {is_data})",
            output=[q.recoil_vis_genboson_p4],
        )
    GenBosonMass = Producer(
        call="lorentzvector::GetMass({df}, {output}, {input})",
        input=[q.recoil_genboson_p4],
        output=[q.genbosonmass],
    )
    with defaults(call=None, input=None, output=None):
        MetBasics = ProducerGroup(
            subproducers=[
                BuildPFMetVector,
                BuildMetVector,
                MetPt_uncorrected,
                MetPhi_uncorrected,
                PFMetPt_uncorrected,
                PFMetPhi_uncorrected,
                MetCov00,
                MetCov01,
                MetCov10,
                MetCov11,
                MetSumEt,
                CalculateGenBosonVector,
                CalculateVisGenBosonVector,
                GenBosonMass,
            ],
        )

with defaults(scopes=["et", "mt", "tt", "em", "mm", "ee"]):
    with defaults(
        call="lorentzvector::PropagateToMET({df}, {output}, {input}, {propagateLeptons})"
    ):
        PropagateLeptonsToMet = Producer(
            input=[q.met_p4, q.p4_1_uncorrected, q.p4_2_uncorrected, q.p4_1, q.p4_2],
            output=[q.met_p4_leptoncorrected],
        )
        PropagateLeptonsToPFMet = Producer(
            input=[q.pfmet_p4, q.p4_1_uncorrected, q.p4_2_uncorrected, q.p4_1, q.p4_2],
            output=[q.pfmet_p4_leptoncorrected],
        )
    with defaults(
        call="physicsobject::PropagateToMET({df}, {output}, {input}, {propagateJets}, {min_jetpt_met_propagation})"
    ):
        PropagateJetsToMet = Producer(
            input=[
                q.met_p4_leptoncorrected,
                q.Jet_pt_corrected,
                nanoAOD.Jet_eta,
                nanoAOD.Jet_phi,
                q.Jet_mass_corrected,
                nanoAOD.Jet_pt,
                nanoAOD.Jet_eta,
                nanoAOD.Jet_phi,
                nanoAOD.Jet_mass,
            ],
            output=[q.met_p4_jetcorrected],
        )
        PropagateJetsToPFMet = Producer(
            input=[
                q.pfmet_p4_leptoncorrected,
                q.Jet_pt_corrected,
                nanoAOD.Jet_eta,
                nanoAOD.Jet_phi,
                q.Jet_mass_corrected,
                nanoAOD.Jet_pt,
                nanoAOD.Jet_eta,
                nanoAOD.Jet_phi,
                nanoAOD.Jet_mass,
            ],
            output=[q.pfmet_p4_jetcorrected],
        )
    with defaults(
        call="""met::RecoilCorrection(
            {df},
            {output},
            {input},
            "{recoil_corrections_file}",
            "{recoil_systematics_file}",
            {applyRecoilCorrections},
            {apply_recoil_resolution_systematic},
            {apply_recoil_response_systematic},
            {recoil_systematic_shift_up},
            {recoil_systematic_shift_down},
            {is_wjets})
        """
    ):
        ApplyRecoilCorrections = Producer(
            input=[
                q.met_p4_jetcorrected,
                q.recoil_genboson_p4,
                q.recoil_vis_genboson_p4,
                q.Jet_pt_corrected,
            ],
            output=[q.met_p4_recoilcorrected],
        )
        ApplyRecoilCorrectionsPFMet = Producer(
            input=[
                q.pfmet_p4_jetcorrected,
                q.recoil_genboson_p4,
                q.recoil_vis_genboson_p4,
                q.Jet_pt_corrected,
            ],
            output=[q.pfmet_p4_recoilcorrected],
        )
    with defaults(call="lorentzvector::GetPt({df}, {output}, {input})"):
        MetPt = Producer(
            input=[q.met_p4_recoilcorrected],
            output=[q.met],
        )
        PFMetPt = Producer(
            input=[q.pfmet_p4_recoilcorrected],
            output=[q.pfmet],
        )
    with defaults(call="lorentzvector::GetPhi({df}, {output}, {input})"):
        MetPhi = Producer(
            input=[q.met_p4_recoilcorrected],
            output=[q.metphi],
        )
        PFMetPhi = Producer(
            input=[q.pfmet_p4_recoilcorrected],
            output=[q.pfmetphi],
        )
    with defaults(call=None, input=None, output=None):
        MetCorrections = ProducerGroup(
            subproducers=[
                PropagateLeptonsToMet,
                PropagateJetsToMet,
                ApplyRecoilCorrections,
                MetPt,
                MetPhi,
            ],
        )
        PFMetCorrections = ProducerGroup(
            subproducers=[
                PropagateLeptonsToPFMet,
                PropagateJetsToPFMet,
                ApplyRecoilCorrectionsPFMet,
                PFMetPt,
                PFMetPhi,
            ],
        )
