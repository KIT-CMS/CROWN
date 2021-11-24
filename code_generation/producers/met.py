import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers used for contruction of met related quantities
####################

BuildMetVector = Producer(
    name="BuildMetVector",
    call="lorentzvectors::buildMet({df}, {input}, {output})",
    input=[
        nanoAOD.MET_pt,
        nanoAOD.MET_phi,
    ],
    output=[q.met_p4],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
BuildPFMetVector = Producer(
    name="BuildPFMetVector",
    call="lorentzvectors::buildMet({df}, {input}, {output})",
    input=[
        nanoAOD.PFMET_pt,
        nanoAOD.PFMET_phi,
    ],
    output=[q.pfmet_p4],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetCov00 = Producer(
    name="MetCov00",
    call="basefunctions::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_covXX,
    ],
    output=[q.metcov00],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetCov01 = Producer(
    name="MetCov01",
    call="basefunctions::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_covXY,
    ],
    output=[q.metcov01],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetCov10 = Producer(
    name="MetCov10",
    call="basefunctions::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_covXY,
    ],
    output=[q.metcov10],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetCov11 = Producer(
    name="MetCov11",
    call="basefunctions::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_covYY,
    ],
    output=[q.metcov11],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetSumEt = Producer(
    name="MetSumEt",
    call="basefunctions::rename<float>({df}, {input}, {output})",
    input=[
        nanoAOD.MET_sumEt,
    ],
    output=[q.metSumEt],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
PropagateLeptonsToMet = Producer(
    name="PropagateLeptonsToMet",
    call="met::propagateLeptonsToMet({df}, {input}, {output}, {propagateLeptons})",
    input=[q.met_p4, q.p4_1_uncorrected, q.p4_2_uncorrected, q.p4_1, q.p4_2],
    output=[q.met_p4_leptoncorrected],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
PropagateJetsToMet = Producer(
    name="PropagateJetsToMet",
    call="met::propagateJetsToMet({df}, {input}, {output}, {propagateJets}, {min_jetpt_met_propagation})",
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
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
CalculateGenBosonVector = Producer(
    name="calculateGenBosonVector",
    call="met::calculateGenBosonVector({df}, {input}, {output}, {is_data})",
    input=[
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
        nanoAOD.GenParticle_pdgId,
        nanoAOD.GenParticle_status,
        nanoAOD.GenParticle_statusFlags,
    ],
    output=[q.recoil_genboson_p4],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
ApplyRecoilCorrections = Producer(
    name="ApplyRecoilCorrections",
    call='met::applyRecoilCorrections({df}, {input}, {output}, "{recoil_corrections_file}", "{recoil_systematics_file}", {applyRecoilCorrections}, {apply_recoil_resolution_systematic}, {apply_recoil_response_systematic}, {recoil_systematic_shift_up}, {recoil_systematic_shift_down}, {is_wj})',
    input=[
        q.met_p4_jetcorrected,
        q.recoil_genboson_p4,
        q.Jet_pt_corrected,
    ],
    output=[q.met_p4_recoilcorrected],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetPt = Producer(
    name="MetPt",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.met_p4_recoilcorrected],
    output=[q.met],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetPhi = Producer(
    name="MetPhi",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.met_p4_recoilcorrected],
    output=[q.metphi],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetPt_uncorrected = Producer(
    name="MetPt_uncorrected",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.met_p4],
    output=[q.met_uncorrected],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetPhi_uncorrected = Producer(
    name="MetPhi_uncorrected",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.met_p4],
    output=[q.metphi_uncorrected],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
PFMetPt_uncorrected = Producer(
    name="PFMetPt_uncorrected",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.pfmet_p4],
    output=[q.pfmet_uncorrected],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
PFMetPhi_uncorrected = Producer(
    name="PFMetPhi_uncorrected",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.pfmet_p4],
    output=[q.pfmetphi_uncorrected],
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
)
MetCorrections = ProducerGroup(
    name="MetCorrections",
    call=None,
    input=None,
    output=None,
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
    subproducers=[
        MetCov00,
        MetCov01,
        MetCov10,
        MetCov11,
        MetSumEt,
        CalculateGenBosonVector,
        PropagateLeptonsToMet,
        PropagateJetsToMet,
        ApplyRecoilCorrections,
        MetPt,
        MetPhi,
    ],
)
UncorrectedMet = ProducerGroup(
    name="MetCorrections",
    call=None,
    input=None,
    output=None,
    scopes=["et", "mt", "tt", "em", "mm", "ee"],
    subproducers=[
        BuildPFMetVector,
        BuildMetVector,
        MetPt_uncorrected,
        MetPhi_uncorrected,
        PFMetPt_uncorrected,
        PFMetPhi_uncorrected,
    ],
)
