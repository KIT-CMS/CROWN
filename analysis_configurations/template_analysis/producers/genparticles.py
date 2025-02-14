from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup

####################
# Set of producers to get the genParticles from the ditaupair
####################

MMGenPair = Producer(
    name="MMGenPair",
    call="ditau_pairselection::buildgenpair({df}, {input}, {output})",
    input=[q.dileptonpair, nanoAOD.Muon_genPartIdx, nanoAOD.Muon_genPartIdx],
    output=[q.gen_dileptonpair],
    scopes=["mm"],
)

####################
# Set of general producers for Gen DiTauPair Quantities
####################

LVGenParticle1 = Producer(
    name="LVGenParticle1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.gen_dileptonpair,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.gen_p4_1],
    scopes=["mm"],
)
LVGenParticle2 = Producer(
    name="LVGenParticle2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.gen_dileptonpair,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.gen_p4_2],
    scopes=["mm"],
)

gen_pt_1 = Producer(
    name="gen_pt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_pt_1],
    scopes=["mm"],
)
gen_pt_2 = Producer(
    name="gen_pt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_pt_2],
    scopes=["mm"],
)
gen_eta_1 = Producer(
    name="gen_eta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_eta_1],
    scopes=["mm"],
)
gen_eta_2 = Producer(
    name="gen_eta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_eta_2],
    scopes=["mm"],
)
gen_phi_1 = Producer(
    name="gen_phi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_phi_1],
    scopes=["mm"],
)
gen_phi_2 = Producer(
    name="gen_phi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_phi_2],
    scopes=["mm"],
)
gen_mass_1 = Producer(
    name="gen_mass_1",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_mass_1],
    scopes=["mm"],
)
gen_mass_2 = Producer(
    name="gen_mass_2",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_mass_2],
    scopes=["mm"],
)
gen_pdgid_1 = Producer(
    name="gen_pdgid_1",
    call="quantities::pdgid({df}, {output}, 0, {input})",
    input=[q.gen_dileptonpair, nanoAOD.GenParticle_pdgId],
    output=[q.gen_pdgid_1],
    scopes=["mm"],
)
gen_pdgid_2 = Producer(
    name="gen_pdgid_2",
    call="quantities::pdgid({df}, {output}, 1, {input})",
    input=[q.gen_dileptonpair, nanoAOD.GenParticle_pdgId],
    output=[q.gen_pdgid_2],
    scopes=["mm"],
)

UnrollGenMuLV1 = ProducerGroup(
    name="UnrollGenMuLV1",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[gen_pt_1, gen_eta_1, gen_phi_1, gen_mass_1, gen_pdgid_1],
)
UnrollGenMuLV2 = ProducerGroup(
    name="UnrollGenMuLV2",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[gen_pt_2, gen_eta_2, gen_phi_2, gen_mass_2, gen_pdgid_2],
)

gen_mm_pair_mass = Producer(
    name="gen_mm_pair_mass",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.gen_p4_1, q.gen_p4_2],
    output=[q.gen_mm_pair_mass],
    scopes=["mm"],
)

MMGenDiTauPairQuantities = ProducerGroup(
    name="GenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mm"],
    subproducers=[
        MMGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenMuLV1,
        UnrollGenMuLV2,
        gen_mm_pair_mass,
    ],
)
