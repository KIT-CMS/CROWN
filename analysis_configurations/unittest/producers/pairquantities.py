from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup, ExtendedVectorProducer

####################
# Set of general producers for DiTauPair Quantities
####################

pt_1 = Producer(
    name="pt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.pt_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
pt_2 = Producer(
    name="pt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.pt_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
eta_1 = Producer(
    name="eta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.eta_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
eta_2 = Producer(
    name="eta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.eta_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
phi_1 = Producer(
    name="phi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.phi_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
phi_2 = Producer(
    name="phi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.phi_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mass_1 = Producer(
    name="mass_1",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.mass_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mass_2 = Producer(
    name="mass_2",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.mass_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
m_vis = Producer(
    name="m_vis",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.m_vis],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
pt_vis = Producer(
    name="pt_vis",
    call="quantities::pt_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.pt_vis],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
####################
# Set of channel specific producers
####################
muon_dxy_1 = Producer(
    name="muon_dxy_1",
    call="quantities::dxy({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_dxy],
    output=[q.dxy_1],
    scopes=["mt", "mm"],
)
muon_dxy_2 = Producer(
    name="muon_dxy_2",
    call="quantities::dxy({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_dxy],
    output=[q.dxy_2],
    scopes=["em", "mm"],
)
electron_dxy_1 = Producer(
    name="electron_dxy_1",
    call="quantities::dxy({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_dxy],
    output=[q.dxy_1],
    scopes=["et", "ee", "em"],
)
electron_dxy_2 = Producer(
    name="electron_dxy_2",
    call="quantities::dxy({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_dxy],
    output=[q.dxy_2],
    scopes=["ee"],
)
tau_dxy_1 = Producer(
    name="tau_dxy_1",
    call="quantities::dxy({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_dxy],
    output=[q.dxy_1],
    scopes=["tt"],
)
tau_dxy_2 = Producer(
    name="tau_dxy_2",
    call="quantities::dxy({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_dxy],
    output=[q.dxy_2],
    scopes=["mt", "et", "tt"],
)
muon_dz_1 = Producer(
    name="muon_dz_1",
    call="quantities::dz({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_dz],
    output=[q.dz_1],
    scopes=["mt", "mm"],
)
muon_dz_2 = Producer(
    name="muon_dz_2",
    call="quantities::dz({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_dz],
    output=[q.dz_2],
    scopes=["em", "mm"],
)
electron_dz_1 = Producer(
    name="electron_dz_1",
    call="quantities::dz({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_dz],
    output=[q.dz_1],
    scopes=["et", "ee", "em"],
)
electron_dz_2 = Producer(
    name="electron_dz_2",
    call="quantities::dz({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_dz],
    output=[q.dz_2],
    scopes=["ee"],
)
tau_dz_1 = Producer(
    name="tau_dz_1",
    call="quantities::dz({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_dz],
    output=[q.dz_1],
    scopes=["tt"],
)
tau_dz_2 = Producer(
    name="tau_dz_2",
    call="quantities::dz({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_dz],
    output=[q.dz_2],
    scopes=["mt", "et", "tt"],
)
muon_q_1 = Producer(
    name="muon_q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_charge],
    output=[q.q_1],
    scopes=["mt", "mm"],
)
muon_q_2 = Producer(
    name="muon_q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_charge],
    output=[q.q_2],
    scopes=["em", "mm"],
)
electron_q_1 = Producer(
    name="electron_q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_charge],
    output=[q.q_1],
    scopes=["et", "ee", "em"],
)
electron_q_2 = Producer(
    name="electron_q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_charge],
    output=[q.q_2],
    scopes=["ee"],
)
tau_q_1 = Producer(
    name="tau_q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_charge],
    output=[q.q_1],
    scopes=["tt"],
)
tau_q_2 = Producer(
    name="tau_q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_charge],
    output=[q.q_2],
    scopes=["mt", "et", "tt"],
)
muon_iso_1 = Producer(
    name="muon_iso_1",
    call="quantities::isolation({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_iso],
    output=[q.iso_1],
    scopes=["mt", "mm"],
)
muon_iso_2 = Producer(
    name="muon_iso_2",
    call="quantities::isolation({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Muon_iso],
    output=[q.iso_2],
    scopes=["em", "mm"],
)
electron_iso_1 = Producer(
    name="electron_iso_1",
    call="quantities::isolation({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_iso],
    output=[q.iso_1],
    scopes=["et", "ee", "em"],
)
electron_iso_2 = Producer(
    name="electron_iso_2",
    call="quantities::isolation({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Electron_iso],
    output=[q.iso_2],
    scopes=["ee"],
)
tau_iso_1 = Producer(
    name="tau_iso_1",
    call="quantities::isolation({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_IDraw],
    output=[q.iso_1],
    scopes=["tt"],
)
tau_iso_2 = Producer(
    name="tau_iso_2",
    call="quantities::isolation({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_IDraw],
    output=[q.iso_2],
    scopes=["mt", "et", "tt"],
)
tau_decaymode_1 = Producer(
    name="decaymode_1",
    call="quantities::tau::decaymode({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_decayMode],
    output=[q.tau_decaymode_1],
    scopes=["tt"],
)
tau_decaymode_1_notau = Producer(
    name="tau_decaymode_1_notau",
    call="basefunctions::DefineQuantity({df}, {output}, -1)",
    input=[],
    output=[q.tau_decaymode_1],
    scopes=["et", "mt", "em", "ee", "mm"],
)
tau_gen_match_1 = Producer(
    name="gen_match_1",
    call="quantities::tau::genmatch({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_genMatch],
    output=[q.gen_match_1],
    scopes=["tt"],
)
taujet_pt_1 = Producer(
    name="taujet_pt_1",
    call="quantities::tau::matching_jet_pt({df}, {output}, 0, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_associatedJet, nanoAOD.Jet_pt],
    output=[q.taujet_pt_1],
    scopes=["tt"],
)
VsJetTauIDFlag_1 = ExtendedVectorProducer(
    name="VsJetTauIDFlag_1",
    call="quantities::tau::TauIDFlag({df}, {output}, 0, {input}, {vsjet_tau_id_WPbit})",
    input=[q.dileptonpair, nanoAOD.Tau_ID_vsJet],
    output="tau_1_vsjet_id_outputname",
    scope=["tt"],
    vec_config="vsjet_tau_id",
)
VsEleTauIDFlag_1 = ExtendedVectorProducer(
    name="VsEleTauIDFlag_1",
    call="quantities::tau::TauIDFlag({df}, {output}, 0, {input}, {vsele_tau_id_WPbit})",
    input=[q.dileptonpair, nanoAOD.Tau_ID_vsEle],
    output="tau_1_vsele_id_outputname",
    scope=["tt"],
    vec_config="vsele_tau_id",
)
VsMuTauIDFlag_1 = ExtendedVectorProducer(
    name="VsMuTauIDFlag_1",
    call="quantities::tau::TauIDFlag({df}, {output}, 0, {input}, {vsmu_tau_id_WPbit})",
    input=[q.dileptonpair, nanoAOD.Tau_ID_vsMu],
    output="tau_1_vsmu_id_outputname",
    scope=["tt"],
    vec_config="vsmu_tau_id",
)

tau_decaymode_2 = Producer(
    name="taudecaymode_2",
    call="quantities::tau::decaymode({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_decayMode],
    output=[q.tau_decaymode_2],
    scopes=["mt", "et", "tt"],
)
tau_decaymode_2_notau = Producer(
    name="tau_decaymode_2_notau",
    call="basefunctions::DefineQuantity({df}, {output}, -1)",
    input=[],
    output=[q.tau_decaymode_2],
    scopes=["em", "ee", "mm"],
)
tau_gen_match_2 = Producer(
    name="taugen_match_2",
    call="quantities::tau::genmatch({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_genMatch],
    output=[q.gen_match_2],
    scopes=["mt", "et", "tt"],
)
taujet_pt_2 = Producer(
    name="taujet_pt_2",
    call="quantities::tau::matching_jet_pt({df}, {output}, 1, {input})",
    input=[q.dileptonpair, nanoAOD.Tau_associatedJet, nanoAOD.Jet_pt],
    output=[q.taujet_pt_2],
    scopes=["mt", "et", "tt"],
)
VsJetTauIDFlag_2 = ExtendedVectorProducer(
    name="VsJetTauIDFlag_2",
    call="quantities::tau::TauIDFlag({df}, {output}, 1, {input}, {vsjet_tau_id_WPbit})",
    input=[q.dileptonpair, nanoAOD.Tau_ID_vsJet],
    output="tau_2_vsjet_id_outputname",
    scope=["et", "mt", "tt"],
    vec_config="vsjet_tau_id",
)
VsEleTauIDFlag_2 = ExtendedVectorProducer(
    name="VsEleTauIDFlag_2",
    call="quantities::tau::TauIDFlag({df}, {output}, 1, {input}, {vsele_tau_id_WPbit})",
    input=[q.dileptonpair, nanoAOD.Tau_ID_vsEle],
    output="tau_2_vsele_id_outputname",
    scope=["et", "mt", "tt"],
    vec_config="vsele_tau_id",
)
VsMuTauIDFlag_2 = ExtendedVectorProducer(
    name="VsMuTauIDFlag_2",
    call="quantities::tau::TauIDFlag({df}, {output}, 1, {input}, {vsmu_tau_id_WPbit})",
    input=[q.dileptonpair, nanoAOD.Tau_ID_vsMu],
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
        tau_gen_match_1,
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
        tau_gen_match_2,
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

Pzetamissvis = Producer(
    name="Pzetamissvis",
    call="quantities::pzetamissvis({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
    output=[q.pzetamissvis],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mTdileptonMET = Producer(
    name="mTdileptonMET",
    call="quantities::mTdileptonMET({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
    output=[q.mTdileptonMET],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mt_1 = Producer(
    name="mt_1",
    call="quantities::mT({df}, {output}, {input})",
    input=[q.p4_1, q.met_p4_recoilcorrected],
    output=[q.mt_1],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mt_2 = Producer(
    name="mt_2",
    call="quantities::mT({df}, {output}, {input})",
    input=[q.p4_2, q.met_p4_recoilcorrected],
    output=[q.mt_2],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
pt_tt = Producer(
    name="pt_tt",
    call="quantities::pt_tt({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
    output=[q.pt_tt],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
pt_ttjj = Producer(
    name="pt_ttjj",
    call="quantities::pt_ttjj({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.jet_p4_1, q.jet_p4_2, q.met_p4_recoilcorrected],
    output=[q.pt_ttjj],
    scopes=["mt", "et", "tt", "em", "ee", "mm"],
)
mt_tot = Producer(
    name="mt_tot",
    call="quantities::mt_tot({df}, {output}, {input})",
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
    subproducers=[Pzetamissvis, mTdileptonMET, mt_1, mt_2, pt_tt, pt_ttjj, mt_tot],
)
p4_fastmtt_mt = Producer(
    name="p4_fastmtt_mt",
    call='quantities::p4_fastmtt({df}, {output}, {input}, "mt")',
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
    call='quantities::p4_fastmtt({df}, {output}, {input}, "et")',
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
    call='quantities::p4_fastmtt({df}, {output}, {input}, "tt")',
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
    call='quantities::p4_fastmtt({df}, {output}, {input}, "em")',
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
    call="quantities::pt({df}, {output}, {input})",
    input=[q.p4_fastmtt],
    output=[q.pt_fastmtt],
    scopes=["mt", "et", "tt", "em"],
)
eta_fastmtt = Producer(
    name="eta_fastmtt",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.p4_fastmtt],
    output=[q.eta_fastmtt],
    scopes=["mt", "et", "tt", "em"],
)
phi_fastmtt = Producer(
    name="phi_fastmtt",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.p4_fastmtt],
    output=[q.phi_fastmtt],
    scopes=["mt", "et", "tt", "em"],
)
m_fastmtt = Producer(
    name="m_fastmtt",
    call="quantities::mass({df}, {output}, {input})",
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
