import code_generation.quantities.output as q
import code_generation.quantities.nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup


####################
# Set of producers to get the genParticles from the ditaupair
####################
MTGenPair = Producer(
    name="GenPair",
    call="pairselection::buildgenpair({df}, {input}, {output})",
    input=[q.ditaupair, nanoAOD.Muon_indexToGen, nanoAOD.Tau_indexToGen],
    output=[q.gen_ditaupair],
    scopes=["mt"],
)
MMGenPair = Producer(
    name="GenPair",
    call="pairselection::buildgenpair({df}, {input}, {output})",
    input=[q.ditaupair, nanoAOD.Muon_indexToGen, nanoAOD.Muon_indexToGen],
    output=[q.gen_ditaupair],
    scopes=["mm"],
)
####################
# Set of general producers for Gen DiTauPair Quantities
####################

LVGenParticle1 = Producer(
    name="LVGenParticle1",
    call="lorentzvectors::build({df}, {input_vec}, 0, {output})",
    input=[
        q.gen_ditaupair,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.gen_p4_1],
    scopes=["mt", "et", "tt", "em", "mm"],
)
LVGenParticle2 = Producer(
    name="LVGenParticle2",
    call="lorentzvectors::build({df}, {input_vec}, 1, {output})",
    input=[
        q.gen_ditaupair,
        nanoAOD.GenParticle_pt,
        nanoAOD.GenParticle_eta,
        nanoAOD.GenParticle_phi,
        nanoAOD.GenParticle_mass,
    ],
    output=[q.gen_p4_2],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_pt_1 = Producer(
    name="gen_pt_1",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_pt_1],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_pt_2 = Producer(
    name="gen_pt_2",
    call="quantities::pt({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_pt_2],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_eta_1 = Producer(
    name="gen_eta_1",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_eta_1],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_eta_2 = Producer(
    name="gen_eta_2",
    call="quantities::eta({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_eta_2],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_phi_1 = Producer(
    name="gen_phi_1",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_phi_1],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_phi_2 = Producer(
    name="gen_phi_2",
    call="quantities::phi({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_phi_2],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_mass_1 = Producer(
    name="gen_mass_1",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.gen_p4_1],
    output=[q.gen_mass_1],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_mass_2 = Producer(
    name="gen_mass_2",
    call="quantities::mass({df}, {output}, {input})",
    input=[q.gen_p4_2],
    output=[q.gen_mass_2],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_pdgid_1 = Producer(
    name="gen_pdgid_1",
    call="quantities::pdgid({df}, {output}, 0, {input})",
    input=[q.gen_ditaupair, nanoAOD.GenParticle_pdgId],
    output=[q.gen_pdgid_1],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_pdgid_2 = Producer(
    name="gen_pdgid_2",
    call="quantities::pdgid({df}, {output}, 1, {input})",
    input=[q.gen_ditaupair, nanoAOD.GenParticle_pdgId],
    output=[q.gen_pdgid_2],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_m_vis = Producer(
    name="gen_m_vis",
    call="quantities::m_vis({df}, {output}, {input_vec})",
    input=[q.gen_p4_1, q.gen_p4_2],
    output=[q.gen_m_vis],
    scopes=["mt", "et", "tt", "em", "mm"],
)
gen_match_2 = Producer(
    name="gen_match_2",
    call="quantities::tau::genmatch({df}, {output}, 1, {input})",
    input=[q.ditaupair, nanoAOD.Tau_genMatch],
    output=[q.gen_match_2],
    scopes=["mt", "et", "tt"],
)
gen_taujet_pt_1 = Producer(
    name="gen_taujet_pt_1",
    call="quantities::tau::matching_genjet_pt({df}, {output}, 1, {input})",
    input=[
        q.ditaupair,
        nanoAOD.Tau_associatedJet,
        nanoAOD.Jet_associatedGenJet,
        nanoAOD.GenJet_pt,
    ],
    output=[q.gen_taujet_pt_1],
    scopes=["tt"],
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
UnrollGenMuLV1 = ProducerGroup(
    name="UnrollGenMuLV1",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "mm"],
    subproducers=[gen_pt_1, gen_eta_1, gen_phi_1, gen_mass_1, gen_pdgid_1],
)
UnrollGenMuLV2 = ProducerGroup(
    name="UnrollGenMuLV2",
    call=None,
    input=None,
    output=None,
    scopes=["em", "mm"],
    subproducers=[gen_pt_2, gen_eta_2, gen_phi_2, gen_mass_2, gen_pdgid_2],
)
UnrollGenElLV1 = ProducerGroup(
    name="UnrollGenElLV1",
    call=None,
    input=None,
    output=None,
    scopes=["em", "ee" "et"],
    subproducers=[gen_pt_1, gen_eta_1, gen_phi_1, gen_mass_1, gen_pdgid_1],
)
UnrollGenElLV2 = ProducerGroup(
    name="UnrollGenElLV2",
    call=None,
    input=None,
    output=None,
    scopes=["ee"],
    subproducers=[gen_pt_2, gen_eta_2, gen_phi_2, gen_mass_2, gen_pdgid_2],
)
UnrollGenTauLV1 = ProducerGroup(
    name="UnrollGenTauLV1",
    call=None,
    input=None,
    output=None,
    scopes=["tt"],
    subproducers=[
        gen_pt_1,
        gen_eta_1,
        gen_phi_1,
        gen_mass_1,
        gen_pdgid_1,
        gen_taujet_pt_1,
    ],
)
UnrollGenTauLV2 = ProducerGroup(
    name="UnrollGenLV2",
    call=None,
    input=None,
    output=None,
    scopes=["mt", "et", "tt"],
    subproducers=[
        gen_pt_2,
        gen_eta_2,
        gen_phi_2,
        gen_mass_2,
        gen_pdgid_2,
        gen_taujet_pt_2,
    ],
)

MTGenDiTauPairQuantities = ProducerGroup(
    name="MTGenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["mt"],
    subproducers=[
        MTGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenMuLV1,
        UnrollGenTauLV2,
        gen_m_vis,
    ],
)
ETGenDiTauPairQuantities = ProducerGroup(
    name="ETGenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["et"],
    subproducers=[
        MTGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenElLV1,
        UnrollGenTauLV2,
        gen_m_vis,
    ],
)
TTGenDiTauPairQuantities = ProducerGroup(
    name="TTGenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["tt"],
    subproducers=[
        MTGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenTauLV1,
        UnrollGenTauLV2,
        gen_m_vis,
    ],
)
EMGenDiTauPairQuantities = ProducerGroup(
    name="EMGenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["em"],
    subproducers=[
        MTGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenElLV1,
        UnrollGenMuLV2,
        gen_m_vis,
    ],
)
EEGenDiTauPairQuantities = ProducerGroup(
    name="EEGenDiTauPairQuantities",
    call=None,
    input=None,
    output=None,
    scopes=["ee"],
    subproducers=[
        MTGenPair,
        LVGenParticle1,
        LVGenParticle2,
        UnrollGenElLV1,
        UnrollGenElLV2,
        gen_m_vis,
    ],
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
        gen_m_vis,
    ],
)
