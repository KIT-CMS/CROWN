from ..quantities import output as q
from ..quantities import nanoAOD as nanoAOD
from code_generation.producer import Producer, Filter
from code_generation.helpers import defaults

####################
# Set of producers used for contruction of MT good pairs and the coressponding lorentz vectors
####################

_nanoAOD_kinematic_muon = [nanoAOD.Muon_pt, nanoAOD.Muon_eta, nanoAOD.Muon_phi, nanoAOD.Muon_mass]

with defaults(scopes=["mm"]):
    with defaults(input=_nanoAOD_kinematic_muon + [q.good_muons_mask], output=[q.dileptonpair]):
        # mm pair with highest pt
        MMPairSelection = Producer(
            call="ditau_pairselection::mumu::PairSelection({df}, {input_vec}, {output}, {mm_pair_min_deltaR})",
        )
        # mm pair closest to Z boson mass
        ZMMPairSelection = Producer(
            call="ditau_pairselection::mumu::ZBosonPairSelection({df}, {input_vec}, {output}, {mm_pair_min_deltaR})",
        )

    GoodMMPairFlag = Producer(
        call="ditau_pairselection::flagGoodPairs({df}, {output}, {input})",
        input=[q.dileptonpair],
        output=[],
    )
    GoodMMPairFilter = Filter(
        call='event::filter::Flags({df}, "GoodMuMuPairs", {input}, "any_of")',
        input=[],
        subproducers=[GoodMMPairFlag],
    )

    LVMu1 = Producer(
        call="lorentzvector::Build({df}, {output}, {input}, 0)",
        input=_nanoAOD_kinematic_muon + [q.dileptonpair],
        output=[q.p4_1],
    )
    LVMu2 = Producer(
        call="lorentzvector::Build({df}, {output}, {input}, 1)",
        input=_nanoAOD_kinematic_muon + [q.dileptonpair],
        output=[q.p4_2],
    )

    # Lorentz vectors for friend trees
    with defaults(call="lorentzvector::Build({df}, {output}, {input})"):
        LVMu1_friend = Producer(
            input=[
                q.pt_1,
                q.eta_1,
                q.phi_1,
                q.mass_1,
            ],
            output=[q.p4_1],
        )
        LVMu2_friend = Producer(
            input=[
                q.pt_2,
                q.eta_2,
                q.phi_2,
                q.mass_2,
            ],
            output=[q.p4_2],
        )
    LV_MM_reconstruction = Producer(
        call="lorentzvector::Sum({df}, {output}, {input})",
        input=[q.p4_1, q.p4_2],
        output=[q.p4_mm_pair],
    )
