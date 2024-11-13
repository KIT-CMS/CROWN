from code_generation.configuration import Configuration
from .producers import all_producers as p
from .quantities import output as q
from .quantities import nanoAOD as nanoAOD

def build_config(
    era,
    sample,
    scopes,
    shifts,
    available_sample_types,
    available_eras,
    available_scopes,
):

    configuration = Configuration(
            era,
            sample,
            scopes,
            shifts,
            available_sample_types,
            available_eras,
            available_scopes,
            global_scope="",
        )
    #########################
    # setup the configuration
    #########################


    # muon base selection:
    configuration.add_config_parameters(
        "mm",
        {
            "n_muon_req": 2,
            "C_muon_req": -1,
        },
    )

    configuration.add_producers(
        "mm",
        [
            p.NumMuonCut,
            p.CMuonCut,
            p.MuonInvMass,
        ],
    )
    
    configuration.add_outputs(
        "mm",
        [
            q.Muon_InvMass,
            # nanoAOD.nMuon,
            # nanoAOD.Muon_charge,
            # nanoAOD.Muon_pt,
            # nanoAOD.Muon_eta,
            # nanoAOD.Muon_phi,
            # nanoAOD.Muon_mass,
        ],
    )

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()