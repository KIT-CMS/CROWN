from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup, ExtendedVectorProducer

####################
# Set of general producers for DiTauPair Quantities
####################

pt_1 = Producer(
    name="pt_1",
    call="lorentzvector::GetPt({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.pt_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
pt_2 = Producer(
    name="pt_2",
    call="lorentzvector::GetPt({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.pt_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
eta_1 = Producer(
    name="eta_1",
    call="lorentzvector::GetEta({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.eta_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
eta_2 = Producer(
    name="eta_2",
    call="lorentzvector::GetEta({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.eta_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
phi_1 = Producer(
    name="phi_1",
    call="lorentzvector::GetPhi({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.phi_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
phi_2 = Producer(
    name="phi_2",
    call="lorentzvector::GetPhi({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.phi_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mass_1 = Producer(
    name="mass_1",
    call="lorentzvector::GetMass({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.mass_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mass_2 = Producer(
    name="mass_2",
    call="lorentzvector::GetMass({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.mass_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
m_vis = Producer(
    name="m_vis",
    call="lorentzvector::GetMass({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2],
    output=[q.m_vis],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
pt_vis = Producer(
    name="pt_vis",
    call="lorentzvector::GetPt({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2],
    output=[q.pt_vis],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
####################
# Set of channel specific producers
####################
muon_dxy_1 = Producer(
    name="muon_dxy_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Muon_dxy, q.dileptonpair],
    output=[q.dxy_1],
    scopes=["mt", "mm"],
)
muon_dxy_2 = Producer(
    name="muon_dxy_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Muon_dxy, q.dileptonpair],
    output=[q.dxy_2],
    scopes=["em", "mm"],
)
electron_dxy_1 = Producer(
    name="electron_dxy_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Electron_dxy, q.dileptonpair],
    output=[q.dxy_1],
    scopes=["et", "ee", "em"],
)
electron_dxy_2 = Producer(
    name="electron_dxy_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Electron_dxy, q.dileptonpair],
    output=[q.dxy_2],
    scopes=["ee"],
)
tau_dxy_1 = Producer(
    name="tau_dxy_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Tau_dxy, q.dileptonpair],
    output=[q.dxy_1],
    scopes=["tt"],
)
tau_dxy_2 = Producer(
    name="tau_dxy_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Tau_dxy, q.dileptonpair],
    output=[q.dxy_2],
    scopes=["mt", "et", "tt"],
)
muon_dz_1 = Producer(
    name="muon_dz_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Muon_dz, q.dileptonpair],
    output=[q.dz_1],
    scopes=["mt", "mm"],
)
muon_dz_2 = Producer(
    name="muon_dz_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Muon_dz, q.dileptonpair],
    output=[q.dz_2],
    scopes=["em", "mm"],
)
electron_dz_1 = Producer(
    name="electron_dz_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Electron_dz, q.dileptonpair],
    output=[q.dz_1],
    scopes=["et", "ee", "em"],
)
electron_dz_2 = Producer(
    name="electron_dz_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Electron_dz, q.dileptonpair],
    output=[q.dz_2],
    scopes=["ee"],
)
tau_dz_1 = Producer(
    name="tau_dz_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Tau_dz, q.dileptonpair],
    output=[q.dz_1],
    scopes=["tt"],
)
tau_dz_2 = Producer(
    name="tau_dz_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Tau_dz, q.dileptonpair],
    output=[q.dz_2],
    scopes=["mt", "et", "tt"],
)
muon_q_1 = Producer(
    name="muon_q_1",
    call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Muon_charge, q.dileptonpair],
    output=[q.q_1],
    scopes=["mt", "mm"],
)
muon_q_2 = Producer(
    name="muon_q_2",
    call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Muon_charge, q.dileptonpair],
    output=[q.q_2],
    scopes=["em", "mm"],
)
electron_q_1 = Producer(
    name="electron_q_1",
    call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Electron_charge, q.dileptonpair],
    output=[q.q_1],
    scopes=["et", "ee", "em"],
)
electron_q_2 = Producer(
    name="electron_q_2",
    call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Electron_charge, q.dileptonpair],
    output=[q.q_2],
    scopes=["ee"],
)
tau_q_1 = Producer(
    name="tau_q_1",
    call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Tau_charge, q.dileptonpair],
    output=[q.q_1],
    scopes=["tt"],
)
tau_q_2 = Producer(
    name="tau_q_2",
    call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Tau_charge, q.dileptonpair],
    output=[q.q_2],
    scopes=["mt", "et", "tt"],
)
muon_iso_1 = Producer(
    name="muon_iso_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Muon_iso, q.dileptonpair],
    output=[q.iso_1],
    scopes=["mt", "mm"],
)
muon_iso_2 = Producer(
    name="muon_iso_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Muon_iso, q.dileptonpair],
    output=[q.iso_2],
    scopes=["em", "mm"],
)
electron_iso_1 = Producer(
    name="electron_iso_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Electron_iso, q.dileptonpair],
    output=[q.iso_1],
    scopes=["et", "ee", "em"],
)
electron_iso_2 = Producer(
    name="electron_iso_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Electron_iso, q.dileptonpair],
    output=[q.iso_2],
    scopes=["ee"],
)
tau_iso_1 = Producer(
    name="tau_iso_1",
    call="event::quantity::Get<float>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Tau_IDraw, q.dileptonpair],
    output=[q.iso_1],
    scopes=["tt"],
)
tau_iso_2 = Producer(
    name="tau_iso_2",
    call="event::quantity::Get<float>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Tau_IDraw, q.dileptonpair],
    output=[q.iso_2],
    scopes=["mt", "et", "tt"],
)
tau_decaymode_1 = Producer(
    name="decaymode_1",
    call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
    input=[nanoAOD.Tau_decayMode, q.dileptonpair],
    output=[q.tau_decaymode_1],
    scopes=["tt"],
)
tau_decaymode_1_notau = Producer(
    name="tau_decaymode_1_notau",
    call="event::quantity::Define({df}, {output}, -1)",
    input=[],
    output=[q.tau_decaymode_1],
    scopes=["et", "mt", "em", "ee", "mm"],
)
taujet_pt_1 = Producer(
    name="taujet_pt_1",
    call="quantities::JetMatching({df}, {output}, {input}, 0)",
    input=[nanoAOD.Jet_pt, nanoAOD.Tau_associatedJet, q.dileptonpair],
    output=[q.taujet_pt_1],
    scopes=["tt"],
)
VsJetTauIDFlag_1 = ExtendedVectorProducer(
    name="VsJetTauIDFlag_1",
    call="physicsobject::tau::quantity::IDFlag({df}, {output}, {input}, 0, {vsjet_tau_id_WPbit})",
    input=[nanoAOD.Tau_ID_vsJet, q.dileptonpair],
    output="tau_1_vsjet_id_outputname",
    scope=["tt"],
    vec_config="vsjet_tau_id",
)
VsEleTauIDFlag_1 = ExtendedVectorProducer(
    name="VsEleTauIDFlag_1",
    call="physicsobject::tau::quantity::IDFlag({df}, {output}, {input}, 0, {vsele_tau_id_WPbit})",
    input=[nanoAOD.Tau_ID_vsEle, q.dileptonpair],
    output="tau_1_vsele_id_outputname",
    scope=["tt"],
    vec_config="vsele_tau_id",
)
VsMuTauIDFlag_1 = ExtendedVectorProducer(
    name="VsMuTauIDFlag_1",
    call="physicsobject::tau::quantity::IDFlag({df}, {output}, {input}, 0, {vsmu_tau_id_WPbit})",
    input=[nanoAOD.Tau_ID_vsMu, q.dileptonpair],
    output="tau_1_vsmu_id_outputname",
    scope=["tt"],
    vec_config="vsmu_tau_id",
)

tau_decaymode_2 = Producer(
    name="taudecaymode_2",
    call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
    input=[nanoAOD.Tau_decayMode, q.dileptonpair],
    output=[q.tau_decaymode_2],
    scopes=["mt", "et", "tt"],
)
tau_decaymode_2_notau = Producer(
    name="tau_decaymode_2_notau",
    call="event::quantity::Define({df}, {output}, -1)",
    input=[],
    output=[q.tau_decaymode_2],
    scopes=["em", "ee", "mm"],
)
taujet_pt_2 = Producer(
    name="taujet_pt_2",
    call="quantities::JetMatching({df}, {output}, {input}, 1)",
    input=[nanoAOD.Jet_pt, nanoAOD.Tau_associatedJet, q.dileptonpair],
    output=[q.taujet_pt_2],
    scopes=["mt", "et", "tt"],
)
VsJetTauIDFlag_2 = ExtendedVectorProducer(
    name="VsJetTauIDFlag_2",
    call="physicsobject::tau::quantity::IDFlag({df}, {output}, {input}, 1, {vsjet_tau_id_WPbit})",
    input=[nanoAOD.Tau_ID_vsJet, q.dileptonpair],
    output="tau_2_vsjet_id_outputname",
    scope=["et", "mt", "tt"],
    vec_config="vsjet_tau_id",
)
VsEleTauIDFlag_2 = ExtendedVectorProducer(
    name="VsEleTauIDFlag_2",
    call="physicsobject::tau::quantity::IDFlag({df}, {output}, {input}, 1, {vsele_tau_id_WPbit})",
    input=[nanoAOD.Tau_ID_vsEle, q.dileptonpair],
    output="tau_2_vsele_id_outputname",
    scope=["et", "mt", "tt"],
    vec_config="vsele_tau_id",
)
VsMuTauIDFlag_2 = ExtendedVectorProducer(
    name="VsMuTauIDFlag_2",
    call="physicsobject::tau::quantity::IDFlag({df}, {output}, {input}, 1, {vsmu_tau_id_WPbit})",
    input=[nanoAOD.Tau_ID_vsMu, q.dileptonpair],
    output="tau_2_vsmu_id_outputname",
    scope=["et", "mt", "tt"],
    vec_config="vsmu_tau_id",
)

UnrollMuLV1 = ProducerGroup(
    name="UnrollMuLV1",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "mm"],
    subproducers=[
        pt_1,
        eta_1,
        phi_1,
        mass_1,
        muon_dxy_1,
        muon_dz_1,
        muon_q_1,
        muon_iso_1,
    ],
)
UnrollMuLV2 = ProducerGroup(
    name="UnrollMuLV2",
    call=None,
    input=None,
    output=None,
    scopes=["mm", "em"],
    subproducers=[
        pt_2,
        eta_2,
        phi_2,
        mass_2,
        muon_dxy_2,
        muon_dz_2,
        muon_q_2,
        muon_iso_2,
    ],
)
UnrollElLV1 = ProducerGroup(
    name="UnrollElLV1",
    call=None,
    input=None,
    output=None,
    scopes=["et", "ee", "em"],
    subproducers=[
        pt_1,
        eta_1,
        phi_1,
        mass_1,
        electron_dxy_1,
        electron_dz_1,
        electron_q_1,
        electron_iso_1,
    ],
)
UnrollElLV2 = ProducerGroup(
    name="UnrollElLV2",
    call=None,
    input=None,
    output=None,
    scopes=["ee"],
    subproducers=[
        pt_2,
        eta_2,
        phi_2,
        mass_2,
        electron_dxy_2,
        electron_dz_2,
        electron_q_2,
        electron_iso_2,
    ],
)
UnrollTauLV1 = ProducerGroup(
    name="UnrollTauLV1",
    call=None,
    input=None,
    output=None,
    scopes=["tt"],
    subproducers=[
        pt_1,
        eta_1,
        phi_1,
        mass_1,
        tau_dxy_1,
        tau_dz_1,
        tau_q_1,
        tau_iso_1,
        tau_decaymode_1,
        taujet_pt_1,
        VsJetTauIDFlag_1,
        VsEleTauIDFlag_1,
        VsMuTauIDFlag_1,
    ],
)
UnrollTauLV2 = ProducerGroup(
    name="UnrollLV2",
    call=None,
    input=None,
    output=None,
    scopes=["et", "mt", "tt"],
    subproducers=[
        pt_2,
        eta_2,
        phi_2,
        mass_2,
        tau_dxy_2,
        tau_dz_2,
        tau_q_2,
        tau_iso_2,
        tau_decaymode_2,
        taujet_pt_2,
        VsJetTauIDFlag_2,
        VsEleTauIDFlag_2,
        VsMuTauIDFlag_2,
    ],
)
#####################
# Producer Groups
#####################

MTDiTauPairQuantities = ProducerGroup(
    name="MTDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt"],
    subproducers=[UnrollMuLV1, UnrollTauLV2, tau_decaymode_1_notau, m_vis, pt_vis],
)
MMDiTauPairQuantities = ProducerGroup(
    name="MMDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[
        UnrollMuLV1,
        UnrollMuLV2,
        tau_decaymode_1_notau,
        tau_decaymode_2_notau,
        m_vis,
        pt_vis,
    ],
)
ETDiTauPairQuantities = ProducerGroup(
    name="ETDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["et"],
    subproducers=[UnrollElLV1, UnrollTauLV2, tau_decaymode_1_notau, m_vis, pt_vis],
)
TTDiTauPairQuantities = ProducerGroup(
    name="TTDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["tt"],
    subproducers=[UnrollTauLV1, UnrollTauLV2, m_vis, pt_vis],
)
EMDiTauPairQuantities = ProducerGroup(
    name="EMDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["em"],
    subproducers=[UnrollElLV1, UnrollMuLV2, m_vis, pt_vis],
)
## advanced event quantities (can be caluculated when ditau pair and met and all jets are determined)
## leptons: q.p4_1, q.p4_2
## met: met_p4_recoilcorrected
## jets: good_jet_collection (if only the leading two are needed: q.jet_p4_1, q.jet_p4_2
## bjets: gen_bjet_collection

LV_dilepton_pair = Producer(
    name="LV_dilepton_pair",
    call="lorentzvector::Sum({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2],
    output=[q.p4_dilepton],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
Pzetamissvis = Producer(
    name="Pzetamissvis",
    call="quantities::PzetaMissVis({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
    output=[q.pzetamissvis],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mTdileptonMET = Producer(
    name="mTdileptonMET",
    call="quantities::TransverseMass({df}, {output}, {input})",
    input=[q.p4_dilepton, q.met_p4_recoilcorrected],
    output=[q.mTdileptonMET],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mt_1 = Producer(
    name="mt_1",
    call="quantities::TransverseMass({df}, {output}, {input})",
    input=[q.p4_1, q.met_p4_recoilcorrected],
    output=[q.mt_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mt_2 = Producer(
    name="mt_2",
    call="quantities::TransverseMass({df}, {output}, {input})",
    input=[q.p4_2, q.met_p4_recoilcorrected],
    output=[q.mt_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
pt_tt = Producer(
    name="pt_tt",
    call="lorentzvector::GetPt({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
    output=[q.pt_tt],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
pt_ttjj = Producer(
    name="pt_ttjj",
    call="lorentzvector::GetPt({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.jet_p4_1, q.jet_p4_2, q.met_p4_recoilcorrected],
    output=[q.pt_ttjj],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mt_tot = Producer(
    name="mt_tot",
    call="quantities::TransverseMass({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
    output=[q.mt_tot],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
DiTauPairMETQuantities = ProducerGroup(
    name="DiTauPairMETQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
    subproducers=[
        LV_dilepton_pair,
        Pzetamissvis,
        mTdileptonMET,
        mt_1,
        mt_2,
        pt_tt,
        pt_ttjj,
        mt_tot,
    ],
)
p4_fastmtt_mt = Producer(
    name="p4_fastmtt_mt",
    call='quantities::FastMtt({df}, {output}, {input}, "mt")',
    input=[
        q.pt_1,
        q.pt_2,
        q.eta_1,
        q.eta_2,
        q.phi_1,
        q.phi_2,
        q.mass_1,
        q.mass_2,
        q.met,
        q.metphi,
        q.metcov00,
        q.metcov01,
        q.metcov11,
        q.tau_decaymode_1,
        q.tau_decaymode_2,
    ],
    output=[q.p4_fastmtt],
    scopes=["mt"],
)
p4_fastmtt_et = Producer(
    name="p4_fastmtt_et",
    call='quantities::FastMtt({df}, {output}, {input}, "et")',
    input=[
        q.pt_1,
        q.pt_2,
        q.eta_1,
        q.eta_2,
        q.phi_1,
        q.phi_2,
        q.mass_1,
        q.mass_2,
        q.met,
        q.metphi,
        q.metcov00,
        q.metcov01,
        q.metcov11,
        q.tau_decaymode_1,
        q.tau_decaymode_2,
    ],
    output=[q.p4_fastmtt],
    scopes=["et"],
)
p4_fastmtt_tt = Producer(
    name="p4_fastmtt_tt",
    call='quantities::FastMtt({df}, {output}, {input}, "tt")',
    input=[
        q.pt_1,
        q.pt_2,
        q.eta_1,
        q.eta_2,
        q.phi_1,
        q.phi_2,
        q.mass_1,
        q.mass_2,
        q.met,
        q.metphi,
        q.metcov00,
        q.metcov01,
        q.metcov11,
        q.tau_decaymode_1,
        q.tau_decaymode_2,
    ],
    output=[q.p4_fastmtt],
    scopes=["tt"],
)
p4_fastmtt_em = Producer(
    name="p4_fastmtt_em",
    call='quantities::FastMtt({df}, {output}, {input}, "em")',
    input=[
        q.pt_1,
        q.pt_2,
        q.eta_1,
        q.eta_2,
        q.phi_1,
        q.phi_2,
        q.mass_1,
        q.mass_2,
        q.met,
        q.metphi,
        q.metcov00,
        q.metcov01,
        q.metcov11,
        q.tau_decaymode_1,
        q.tau_decaymode_2,
    ],
    output=[q.p4_fastmtt],
    scopes=["em"],
)

pt_fastmtt = Producer(
    name="pt_fastmtt",
    call="lorentzvector::GetPt({df}, {output}, {input})",
    input=[q.p4_fastmtt],
    output=[q.pt_fastmtt],
    scopes=["mt", "et", "tt", "em"],
)
eta_fastmtt = Producer(
    name="eta_fastmtt",
    call="lorentzvector::GetEta({df}, {output}, {input})",
    input=[q.p4_fastmtt],
    output=[q.eta_fastmtt],
    scopes=["mt", "et", "tt", "em"],
)
phi_fastmtt = Producer(
    name="phi_fastmtt",
    call="lorentzvector::GetPhi({df}, {output}, {input})",
    input=[q.p4_fastmtt],
    output=[q.phi_fastmtt],
    scopes=["mt", "et", "tt", "em"],
)
m_fastmtt = Producer(
    name="m_fastmtt",
    call="lorentzvector::GetMass({df}, {output}, {input})",
    input=[q.p4_fastmtt],
    output=[q.m_fastmtt],
    scopes=["mt", "et", "tt", "em"],
)
FastMTTQuantities = ProducerGroup(
    name="FastMTTQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "et", "tt", "em"],
    subproducers={
        "mt": [p4_fastmtt_mt, pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt],
        "et": [p4_fastmtt_et, pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt],
        "tt": [p4_fastmtt_tt, pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt],
        "em": [p4_fastmtt_em, pt_fastmtt, eta_fastmtt, phi_fastmtt, m_fastmtt],
    },
)
