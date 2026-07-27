from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, ProducerGroup
from code_generation.helpers import defaults

####################
# Set of general producers for DiTauPair Quantities
####################

with defaults(scopes=["mm"]):
    with defaults(input=[q.p4_1]):
        pt_1 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.pt_1],
        )
        eta_1 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.eta_1],
        )
        phi_1 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.phi_1],
        )
        mass_1 = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.mass_1],
        )
    with defaults(input=[q.p4_2]):
        pt_2 = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.pt_2],
        )
        eta_2 = Producer(
            call="lorentzvector::GetEta({df}, {output}, {input})",
            output=[q.eta_2],
        )
        phi_2 = Producer(
            call="lorentzvector::GetPhi({df}, {output}, {input})",
            output=[q.phi_2],
        )
        mass_2 = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.mass_2],
        )
    with defaults(input=[nanoAOD.Muon_charge, q.dileptonpair]):
        muon_q_1 = Producer(
            call="event::quantity::Get<int>({df}, {output}, {input}, 0)",
            output=[q.q_1],
        )
        muon_q_2 = Producer(
            call="event::quantity::Get<int>({df}, {output}, {input}, 1)",
            output=[q.q_2],
        )
    with defaults(call=None, input=None, output=None):
        UnrollMuLV1 = ProducerGroup(
            subproducers=[
                pt_1,
                eta_1,
                phi_1,
                mass_1,
                muon_q_1,
            ],
        )
        UnrollMuLV2 = ProducerGroup(
            subproducers=[
                pt_2,
                eta_2,
                phi_2,
                mass_2,
                muon_q_2,
            ],
        )

    ### for friends
    with defaults(input=[q.p4_1, q.p4_2]):
        mm_pair_mass = Producer(
            call="lorentzvector::GetMass({df}, {output}, {input})",
            output=[q.mm_pair_mass],
        )
        mm_pair_pt = Producer(
            call="lorentzvector::GetPt({df}, {output}, {input})",
            output=[q.mm_pair_pt],
        )
    with defaults(call=None, input=None, output=None):
        MMDiTauPairQuantities = ProducerGroup(
            subproducers=[UnrollMuLV1, UnrollMuLV2, mm_pair_mass, mm_pair_pt],
        )
