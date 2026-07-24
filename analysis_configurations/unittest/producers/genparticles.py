from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

####################
# Set of producers to get the genParticles from the ditaupair
####################
with defaults(
    call="ditau_pairselection::buildgenpair({df}, {output}, {input})",
    output=[q.gen_dileptonpair],
):
    MTGenPair = Producer(
        input=[q.dileptonpair, nanoAOD.Muon_genPartIdx, nanoAOD.Tau_genPartIdx],
        scopes=["mt"],
    )
    ETGenPair = Producer(
        input=[q.dileptonpair, nanoAOD.Electron_genPartIdx, nanoAOD.Tau_genPartIdx],
        scopes=["et"],
    )
    TTGenPair = Producer(
        input=[q.dileptonpair, nanoAOD.Tau_genPartIdx, nanoAOD.Tau_genPartIdx],
        scopes=["tt"],
    )
    EMGenPair = Producer(
        input=[q.dileptonpair, nanoAOD.Electron_genPartIdx, nanoAOD.Muon_genPartIdx],
        scopes=["em"],
    )
    MMGenPair = Producer(
        input=[q.dileptonpair, nanoAOD.Muon_genPartIdx, nanoAOD.Muon_genPartIdx],
        scopes=["mm"],
    )
    EEGenPair = Producer(
        input=[
            q.dileptonpair,
            nanoAOD.Electron_genPartIdx,
            nanoAOD.Electron_genPartIdx,
        ],
        scopes=["ee"],
    )
MMTrueGenPair = Producer(
    call="""ditau_pairselection::buildtruegenpair(
        {df},
        {output},
        {input},
        {truegen_mother_pdgid}, 
        {truegen_daughter_1_pdgid}, 
        {truegen_daugher_2_pdgid})
        """,
    input=[
        nanoAOD.GenPart_statusFlags,
        nanoAOD.GenPart_status,
        nanoAOD.GenPart_pdgId,
        nanoAOD.GenPart_genPartIdxMother,
        nanoAOD.GenPart_pt,
    ],
    output=[q.truegenpair],
    scopes=["mm"],
)
####################
# Set of general producers for Gen DiTauPair Quantities
####################
_nanoAOD_kinematic_GenPart = [
    nanoAOD.GenPart_pt,
    nanoAOD.GenPart_eta,
    nanoAOD.GenPart_phi,
    nanoAOD.GenPart_mass,
]
with defaults(scopes=["mt", "et", "tt", "em", "mm", "ee"]):
    with defaults(input=_nanoAOD_kinematic_GenPart + [q.gen_dileptonpair]):
        LVGenParticle1 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 0)",
            output=[q.gen_p4_1],
        )
        LVGenParticle2 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 1)",
            output=[q.gen_p4_2],
        )
    with defaults(input=_nanoAOD_kinematic_GenPart + [q.truegenpair]):
        LVTrueGenParticle1 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 0)",
            output=[q.gen_p4_1],
        )
        LVTrueGenParticle2 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 1)",
            output=[q.gen_p4_2],
        )
    with defaults(input=[q.gen_p4_1]):
        gen_pt_1 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.gen_pt_1],
        )
        gen_eta_1 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.gen_eta_1],
        )
        gen_phi_1 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.gen_phi_1],
        )
        gen_mass_1 = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.gen_mass_1],
        )
    with defaults(input=[q.gen_p4_2]):
        gen_pt_2 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.gen_pt_2],
        )
        gen_eta_2 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.gen_eta_2],
        )
        gen_phi_2 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.gen_phi_2],
        )
        gen_mass_2 = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.gen_mass_2],
        )
    gen_pdgid_1 = Producer(
        call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
        input=[nanoAOD.GenPart_pdgId, q.gen_dileptonpair],
        output=[q.gen_pdgid_1],
    )
    gen_pdgid_2 = Producer(
        call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
        input=[nanoAOD.GenPart_pdgId, q.gen_dileptonpair],
        output=[q.gen_pdgid_2],
    )
    gen_m_vis = Producer(
        call="lorentzvector::GetMass({df}, {output}, {input})",
        input=[q.gen_p4_1, q.gen_p4_2],
        output=[q.gen_m_vis],
    )

gen_taujet_pt_1 = Producer(
    call="event::quantity::GetGenJetForObject<float>({df}, {output}, {input}, 0)",
    input=[
        nanoAOD.GenJet_pt,
        nanoAOD.Jet_genJetIdx,
        nanoAOD.Tau_jetIdx,
        q.dileptonpair,
    ],
    output=[q.gen_taujet_pt_1],
    scopes=["tt"],
)
gen_taujet_pt_2 = Producer(
    call="event::quantity::GetGenJetForObject<float>({df}, {output}, {input}, 1)",
    input=[
        nanoAOD.GenJet_pt,
        nanoAOD.Jet_genJetIdx,
        nanoAOD.Tau_jetIdx,
        q.dileptonpair,
    ],
    output=[q.gen_taujet_pt_2],
    scopes=["mt", "et", "tt"],
)
with defaults(call=None, input=None, output=None):
    UnrollGenMuLV1 = ProducerGroup(
        scopes=["mt", "mm"],
        subproducers=[gen_pt_1, gen_eta_1, gen_phi_1, gen_mass_1, gen_pdgid_1],
    )
    UnrollGenMuLV2 = ProducerGroup(
        scopes=["em", "mm"],
        subproducers=[gen_pt_2, gen_eta_2, gen_phi_2, gen_mass_2, gen_pdgid_2],
    )
    UnrollGenElLV1 = ProducerGroup(
        scopes=["em", "ee", "et"],
        subproducers=[gen_pt_1, gen_eta_1, gen_phi_1, gen_mass_1, gen_pdgid_1],
    )
    UnrollGenElLV2 = ProducerGroup(
        scopes=["ee"],
        subproducers=[gen_pt_2, gen_eta_2, gen_phi_2, gen_mass_2, gen_pdgid_2],
    )
    UnrollGenTauLV1 = ProducerGroup(
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
        scopes=["et"],
        subproducers=[
            ETGenPair,
            LVGenParticle1,
            LVGenParticle2,
            UnrollGenElLV1,
            UnrollGenTauLV2,
            gen_m_vis,
        ],
    )
    TTGenDiTauPairQuantities = ProducerGroup(
        scopes=["tt"],
        subproducers=[
            TTGenPair,
            LVGenParticle1,
            LVGenParticle2,
            UnrollGenTauLV1,
            UnrollGenTauLV2,
            gen_m_vis,
        ],
    )
    EMGenDiTauPairQuantities = ProducerGroup(
        scopes=["em"],
        subproducers=[
            EMGenPair,
            LVGenParticle1,
            LVGenParticle2,
            UnrollGenElLV1,
            UnrollGenMuLV2,
            gen_m_vis,
        ],
    )
    EEGenDiTauPairQuantities = ProducerGroup(
        scopes=["ee"],
        subproducers=[
            EEGenPair,
            LVGenParticle1,
            LVGenParticle2,
            UnrollGenElLV1,
            UnrollGenElLV2,
            gen_m_vis,
        ],
    )
    MMGenDiTauPairQuantities = ProducerGroup(
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
    MMTrueGenDiTauPairQuantities = ProducerGroup(
        scopes=["mm"],
        subproducers=[
            MMTrueGenPair,
            LVTrueGenParticle1,
            LVTrueGenParticle2,
            UnrollGenMuLV1,
            UnrollGenMuLV2,
            gen_m_vis,
        ],
    )

#######################
# DiTau Genmatching
#######################

with defaults(scopes=["mt", "et", "tt", "em", "ee", "mm"]):
    HadronicGenTaus = Producer(
        call="genparticles::tau::HadronicGenTaus({df}, {output}, {input})",
        input=[
            nanoAOD.GenPart_pdgId,
            nanoAOD.GenPart_statusFlags,
            nanoAOD.GenPart_genPartIdxMother,
        ],
        output=[q.hadronic_gen_taus],
    )
    with defaults(call="genparticles::tau::GenMatching({df}, {output}, {input})"):
        GenMatchP1 = Producer(
            input=[
                q.hadronic_gen_taus,
                nanoAOD.GenPart_pdgId,
                nanoAOD.GenPart_statusFlags,
                nanoAOD.GenPart_pt,
                nanoAOD.GenPart_eta,
                nanoAOD.GenPart_phi,
                nanoAOD.GenPart_mass,
                q.p4_1,
            ],
            output=[q.gen_match_1],
        )

        GenMatchP2 = Producer(
            input=[
                q.hadronic_gen_taus,
                nanoAOD.GenPart_pdgId,
                nanoAOD.GenPart_statusFlags,
                nanoAOD.GenPart_pt,
                nanoAOD.GenPart_eta,
                nanoAOD.GenPart_phi,
                nanoAOD.GenPart_mass,
                q.p4_2,
            ],
            output=[q.gen_match_2],
        )
    with defaults(call=None, input=None, output=None):
        GenMatching = ProducerGroup(
            subproducers=[
                HadronicGenTaus,
                GenMatchP1,
                GenMatchP2,
            ],
        )
