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
    scopes=["mt", "et", "tt", "em"],
)
pt_2 = Producer(
    name="pt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.pt_2],
    scopes=["mt", "et", "tt", "em"],
)
eta_1 = Producer(
    name="eta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.eta_1],
    scopes=["mt", "et", "tt", "em"],
)
eta_2 = Producer(
    name="eta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.eta_2],
    scopes=["mt", "et", "tt", "em"],
)
phi_1 = Producer(
    name="phi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.phi_1],
    scopes=["mt", "et", "tt", "em"],
)
phi_2 = Producer(
    name="phi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.phi_2],
    scopes=["mt", "et", "tt", "em"],
)
mass_1 = Producer(
    name="mass_1",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.p4_1],
    output=[q.mass_1],
    scopes=["mt", "et", "tt", "em"],
)
mass_2 = Producer(
    name="mass_2",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.p4_2],
    output=[q.mass_2],
    scopes=["mt", "et", "tt", "em"],
)
m_vis = Producer(
    name="m_vis",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.p4_1, q.p4_2],
    output=[q.m_vis],
    scopes=["mt", "et", "tt", "em"],
)
####################
# Set of channel specific producers
####################
dxy_1 = Producer(
    name="dxy_1",
    call="quantities::dxy({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_dxy],
    output=[q.dxy_1],
    scopes=["mt"],
)
dxy_2 = Producer(
    name="dxy_2",
    call="quantities::dxy({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_dxy],
    output=[q.dxy_2],
    scopes=["mt", "et", "tt"],
)
dz_1 = Producer(
    name="dz_1",
    call="quantities::dz({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_dz],
    output=[q.dz_1],
    scopes=["mt"],
)
dz_2 = Producer(
    name="dz_2",
    call="quantities::dz({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_dz],
    output=[q.dz_2],
    scopes=["mt", "et", "tt"],
)
q_1 = Producer(
    name="q_1",
    call="quantities::charge({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_charge],
    output=[q.q_1],
    scopes=["mt"],
)
q_2 = Producer(
    name="q_2",
    call="quantities::charge({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_charge],
    output=[q.q_2],
    scopes=["mt", "et", "tt"],
)
iso_1 = Producer(
    name="iso_1",
    call="quantities::isolation({df}, {output}, 0, {input})",
    input=[q.ditaupair, nanoAOD.Muon_iso],
    output=[q.iso_1],
    scopes=["mt"],
)
iso_2 = Producer(
    name="iso_2",
    call="quantities::isolation({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_IDraw],
    output=[q.iso_2],
    scopes=["mt", "et", "tt"],
)
decaymode_2 = Producer(
    name="decaymode_2",
    call="quantities::tau::decaymode({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_decayMode],
    output=[q.decaymode_2],
    scopes=["mt", "et", "tt"],
)
gen_match_2 = Producer(
    name="gen_match_2",
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
gen_taujet_pt_2 = Producer(
    name="gen_taujet_pt_2",
    call="quantities::tau::matching_genjet_pt({df}, {output}, 1, {input})",
    input=[
        q.ditaupair,
        nanoAOD.Tau_associatedJet,
        nanoAOD.Jet_associatedGenJet,
        nanoAOD.GenJet_pt,
    ],
    output=[q.gen_taujet_pt_2],
    scopes=["mt", "et", "tt"],
)
UnrollLV1 = ProducerGroup(
    name="UnrollLV1",
    call=None,
    input=None,
    output=None,
    scopes=["mt"],
    subproducers=[pt_1, eta_1, phi_1, mass_1, dxy_1, dz_1, q_1, iso_1],
)
UnrollLV2 = ProducerGroup(
    name="UnrollLV2",
    call=None,
    input=None,
    output=None,
    scopes=["mt"],
    subproducers=[
        pt_2,
        eta_2,
        phi_2,
        mass_2,
        dxy_2,
        dz_2,
        q_2,
        iso_2,
        decaymode_2,
        gen_match_2,
        taujet_pt_2,
        gen_taujet_pt_2,
    ],
)
DiTauPairQuantities = ProducerGroup(
    name="DiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "et", "tt", "em"],
    subproducers=[UnrollLV1, UnrollLV2, m_vis],
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
    scopes=["mt", "et", "tt", "em"],
)
mTdileptonMET = Producer(
    name="mTdileptonMET",
    call="quantities::mTdileptonMET({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
    output=[q.mTdileptonMET],
    scopes=["mt", "et", "tt", "em"],
)
mt_1 = Producer(
    name="mt_1",
    call="quantities::mT({df}, {output}, {input})",
    input=[q.p4_1, q.met_p4_recoilcorrected],
    output=[q.mt_1],
    scopes=["mt", "et", "tt", "em"],
)
mt_2 = Producer(
    name="mt_2",
    call="quantities::mT({df}, {output}, {input})",
    input=[q.p4_2, q.met_p4_recoilcorrected],
    output=[q.mt_2],
    scopes=["mt", "et", "tt", "em"],
)
pt_tt = Producer(
    name="pt_tt",
    call="quantities::pt_tt({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.met_p4_recoilcorrected],
    output=[q.pt_tt],
    scopes=["mt", "et", "tt", "em"],
)
pt_ttjj = Producer(
    name="pt_ttjj",
    call="quantities::pt_ttjj({df}, {output}, {input})",
    input=[q.p4_1, q.p4_2, q.jet_p4_1, q.jet_p4_2, q.met_p4_recoilcorrected],
    output=[q.pt_ttjj],
    scopes=["mt", "et", "tt", "em"],
)
DiTauPairMETQuantities = ProducerGroup(
    name="DiTauPairMETQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "et", "tt", "em"],
    subproducers=[Pzetamissvis, mTdileptonMET, mt_1, mt_2, pt_tt, pt_ttjj],
)
