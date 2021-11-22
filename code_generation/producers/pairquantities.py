import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

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
    input=[q.ditaupair, nanoAOD.Muon_dxy],
    output=[q.dxy_1],
    scopes=["mt", "mm"],
)
muon_dxy_2 = Producer(
    name="muon_dxy_2",
    call="quantities::dxy({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Muon_dxy],
    output=[q.dxy_2],
    scopes=["em", "mm"],
)
tau_dxy_1 = Producer(
    name="tau_dxy_1",
    call="quantities::dxy({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Tau_dxy],
    output=[q.dxy_1],
    scopes=["tt"],
)
tau_dxy_2 = Producer(
    name="tau_dxy_2",
    call="quantities::dxy({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_dxy],
    output=[q.dxy_2],
    scopes=["mt", "et", "tt"],
)
muon_dz_1 = Producer(
    name="muon_dz_1",
    call="quantities::dz({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_dz],
    output=[q.dz_1],
    scopes=["mt", "mm"],
)
muon_dz_2 = Producer(
    name="muon_dz_2",
    call="quantities::dz({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Muon_dz],
    output=[q.dz_2],
    scopes=["em", "mm"],
)
tau_dz_1 = Producer(
    name="tau_dz_1",
    call="quantities::dz({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Tau_dz],
    output=[q.dz_1],
    scopes=["tt"],
)
tau_dz_2 = Producer(
    name="tau_dz_2",
    call="quantities::dz({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_dz],
    output=[q.dz_2],
    scopes=["mt", "et", "tt"],
)
muon_q_1 = Producer(
    name="muon_q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_charge],
    output=[q.q_1],
    scopes=["mt", "mm"],
)
muon_q_2 = Producer(
    name="muon_q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Muon_charge],
    output=[q.q_2],
    scopes=["em", "mm"],
)
tau_q_1 = Producer(
    name="tau_q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Tau_charge],
    output=[q.q_1],
    scopes=["tt"],
)
tau_q_2 = Producer(
    name="tau_q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_charge],
    output=[q.q_2],
    scopes=["mt", "et", "tt"],
)
muon_iso_1 = Producer(
    name="muon_iso_1",
    call="quantities::isolation({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_iso],
    output=[q.iso_1],
    scopes=["mt", "mm"],
)
muon_iso_2 = Producer(
    name="muon_iso_2",
    call="quantities::isolation({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Muon_iso],
    output=[q.iso_2],
    scopes=["em", "mm"],
)
tau_iso_1 = Producer(
    name="tau_iso_1",
    call="quantities::isolation({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Tau_IDraw],
    output=[q.iso_1],
    scopes=["tt"],
)
tau_iso_2 = Producer(
    name="tau_iso_2",
    call="quantities::isolation({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_IDraw],
    output=[q.iso_2],
    scopes=["mt", "et", "tt"],
)
tau_decaymode_1 = Producer(
    name="decaymode_1",
    call="quantities::tau::decaymode({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Tau_decayMode],
    output=[q.decaymode_1],
    scopes=["tt"],
)
tau_gen_match_1 = Producer(
    name="gen_match_1",
    call="quantities::tau::genmatch({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Tau_genMatch],
    output=[q.gen_match_1],
    scopes=["tt"],
)
taujet_pt_1 = Producer(
    name="taujet_pt_1",
    call="quantities::tau::matching_jet_pt({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Tau_associatedJet, nanoAOD.Jet_pt],
    output=[q.taujet_pt_1],
    scopes=["tt"],
)

tau_decaymode_2 = Producer(
    name="taudecaymode_2",
    call="quantities::tau::decaymode({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_decayMode],
    output=[q.decaymode_2],
    scopes=["mt", "et", "tt"],
)
tau_gen_match_2 = Producer(
    name="taugen_match_2",
    call="quantities::tau::genmatch({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_genMatch],
    output=[q.gen_match_2],
    scopes=["mt", "et", "tt"],
)
taujet_pt_2 = Producer(
    name="taujet_pt_2",
    call="quantities::tau::matching_jet_pt({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_associatedJet, nanoAOD.Jet_pt],
    output=[q.taujet_pt_2],
    scopes=["mt", "et", "tt"],
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
    ],
)
MTDiTauPairQuantities = ProducerGroup(
    name="DiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt"],
    subproducers=[UnrollMuLV1, UnrollTauLV2, m_vis, pt_vis],
)
MMDiTauPairQuantities = ProducerGroup(
    name="DiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[UnrollMuLV1, UnrollMuLV2, m_vis, pt_vis],
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
