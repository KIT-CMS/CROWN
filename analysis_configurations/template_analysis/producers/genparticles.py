from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

with defaults(scopes=["mm"]):
    MMGenPair = Producer(
        call="ditau_pairselection::buildgenpair({df}, {output}, {input})",
        input=[q.dileptonpair, nanoAOD.Muon_genPartIdx, nanoAOD.Muon_genPartIdx],
        output=[q.gen_dileptonpair],
    )

    with defaults(
        input=[
            nanoAOD.GenPart_pt,
            nanoAOD.GenPart_eta,
            nanoAOD.GenPart_phi,
            nanoAOD.GenPart_mass,
            q.gen_dileptonpair,
        ],
    ):
        LVGenParticle1 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 0)",
            output=[q.gen_p4_1],
        )
        LVGenParticle2 = Producer(
            call="lorentzvector::Build({df}, {output}, {input}, 1)",
            output=[q.gen_p4_2],
        )

    with defaults(input=[q.get_p4_1]):
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

    with defaults(input=[nanoAOD.GenPart_pdgId, q.gen_dileptonpair]):
        gen_pdgid_1 = Producer(
            call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
            output=[q.gen_pdgid_1],
        )
        gen_pdgid_2 = Producer(
            call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
            output=[q.gen_pdgid_2],
        )

    with defaults(input=[q.gen_p4_1, q.gen_p4_2]):
        gen_mm_pair_mass = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.gen_mm_pair_mass],
        )
        LV_GenMM_reconstruction = Producer(
            call="lorentzvector::Sum({df}, {output}, {input})",
            output=[q.gen_p4_mm_pair],
        )

    # Lorentz vectors for friend trees
    with defaults(call="lorentzvector::Build({df}, {output}, {input})"):
        LVGenParticle1_friend = Producer(
            input=[q.gen_pt_1, q.gen_eta_1, q.gen_phi_1, q.gen_mass_1],
            output=[q.gen_p4_1],
        )
        LVGenParticle2_friend = Producer(
            input=[q.gen_pt_2, q.gen_eta_2, q.gen_phi_2, q.gen_mass_2],
            output=[q.gen_p4_2],
        )

    with defaults(call=None, input=None, output=None):
        UnrollGenMuLV1 = ProducerGroup(
            subproducers=[gen_pt_1, gen_eta_1, gen_phi_1, gen_mass_1, gen_pdgid_1],
        )
        UnrollGenMuLV2 = ProducerGroup(
            subproducers=[gen_pt_2, gen_eta_2, gen_phi_2, gen_mass_2, gen_pdgid_2],
        )
        MMGenDiTauPairQuantities = ProducerGroup(
            subproducers=[
                MMGenPair,
                LVGenParticle1,
                LVGenParticle2,
                UnrollGenMuLV1,
                UnrollGenMuLV2,
                gen_mm_pair_mass,
            ],
        )
