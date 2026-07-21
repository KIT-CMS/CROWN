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

    LVGenParticle1 = Producer(
        call="lorentzvector::Build({df}, {output}, {input}, 0)",
        input=[
            nanoAOD.GenPart_pt,
            nanoAOD.GenPart_eta,
            nanoAOD.GenPart_phi,
            nanoAOD.GenPart_mass,
            q.gen_dileptonpair,
        ],
        output=[q.gen_p4_1],
    )
    LVGenParticle2 = Producer(
        call="lorentzvector::Build({df}, {output}, {input}, 1)",
        input=[
            nanoAOD.GenPart_pt,
            nanoAOD.GenPart_eta,
            nanoAOD.GenPart_phi,
            nanoAOD.GenPart_mass,
            q.gen_dileptonpair,
        ],
        output=[q.gen_p4_2],
    )

    gen_pt_1 = Producer(
        call="lorentzvector::GetPt({df}, {output}, {input})",
        input=[q.gen_p4_1],
        output=[q.gen_pt_1],
    )
    gen_pt_2 = Producer(
        call="lorentzvector::GetPt({df}, {output}, {input})",
        input=[q.gen_p4_2],
        output=[q.gen_pt_2],
    )
    gen_eta_1 = Producer(
        call="lorentzvector::GetEta({df}, {output}, {input})",
        input=[q.gen_p4_1],
        output=[q.gen_eta_1],
    )
    gen_eta_2 = Producer(
        call="lorentzvector::GetEta({df}, {output}, {input})",
        input=[q.gen_p4_2],
        output=[q.gen_eta_2],
    )
    gen_phi_1 = Producer(
        call="lorentzvector::GetPhi({df}, {output}, {input})",
        input=[q.gen_p4_1],
        output=[q.gen_phi_1],
    )
    gen_phi_2 = Producer(
        call="lorentzvector::GetPhi({df}, {output}, {input})",
        input=[q.gen_p4_2],
        output=[q.gen_phi_2],
    )
    gen_mass_1 = Producer(
        call="lorentzvector::GetMass({df}, {output}, {input})",
        input=[q.gen_p4_1],
        output=[q.gen_mass_1],
    )
    gen_mass_2 = Producer(
        call="lorentzvector::GetMass({df}, {output}, {input})",
        input=[q.gen_p4_2],
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

    UnrollGenMuLV1 = ProducerGroup(
        call=None,
        input=None,
        output=None,
        subproducers=[gen_pt_1, gen_eta_1, gen_phi_1, gen_mass_1, gen_pdgid_1],
    )
    UnrollGenMuLV2 = ProducerGroup(
        call=None,
        input=None,
        output=None,
        subproducers=[gen_pt_2, gen_eta_2, gen_phi_2, gen_mass_2, gen_pdgid_2],
    )

    gen_mm_pair_mass = Producer(
        call="lorentzvector::GetMass({df}, {output}, {input})",
        input=[q.gen_p4_1, q.gen_p4_2],
        output=[q.gen_mm_pair_mass],
    )

    MMGenDiTauPairQuantities = ProducerGroup(
        call=None,
        input=None,
        output=None,
        subproducers=[
            MMGenPair,
            LVGenParticle1,
            LVGenParticle2,
            UnrollGenMuLV1,
            UnrollGenMuLV2,
            gen_mm_pair_mass,
        ],
    )

    # Lorentz vectors for friend trees
    LVGenParticle1_friend = Producer(
        call="lorentzvector::Build({df}, {output}, {input})",
        input=[
            q.gen_pt_1,
            q.gen_eta_1,
            q.gen_phi_1,
            q.gen_mass_1,
        ],
        output=[q.gen_p4_1],
    )
    LVGenParticle2_friend = Producer(
        call="lorentzvector::Build({df}, {output}, {input})",
        input=[
            q.gen_pt_2,
            q.gen_eta_2,
            q.gen_phi_2,
            q.gen_mass_2,
        ],
        output=[q.gen_p4_2],
    )
    LV_GenMM_reconstruction = Producer(
        call="lorentzvector::Sum({df}, {output}, {input})",
        input=[q.gen_p4_1, q.gen_p4_2],
        output=[q.gen_p4_mm_pair],
    )
