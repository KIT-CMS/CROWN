from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List

from .producers import muon_sf_friends as muon_sf_friends
from .producers import event as event
from .producers import genparticles as genparticles
from .producers import jets as jets
from .producers import met as met
from .producers import muons as muons
from .producers import pairquantities as pairquantities
from .producers import pairselection as pairselection
from .producers import scalefactors as scalefactors
from .producers import taus as taus
from .producers import triggers as triggers
from .quantities import nanoAOD as nanoAOD
from .quantities import output as q
# from code_generation.configuration import Configuration
from code_generation.friend_trees import FriendTreeConfiguration


def build_config(
    era: str,
    sample: str,
    scopes: List[str],
    shifts: List[str],
    available_sample_types: List[str],
    available_eras: List[str],
    available_scopes: List[str],
):
    input_testfile = "/work/sbrommer/ntuple_prototype/CROWNTestingSamples/CrownNtuple.root"
    run_nominal = (
        True if "nominal" in shifts else False
    )
    configuration = FriendTreeConfiguration(
        era,
        sample,
        scopes,
        shifts,
        available_sample_types,
        available_eras,
        available_scopes,
        input_testfile,
        run_nominal=run_nominal,
    )

    configuration.add_config_parameters(
        ["mt"],
        {
            "muon_sf_file": "data/embedding/muon_2018UL.json.gz",
            "muon_id_sf": "ID_pt_eta_bins",
            "muon_iso_sf": "Iso_pt_eta_bins",
        },
    )

    configuration.add_producers(
        ["mt"],
        [muon_sf_friends.MuonIDSF_friends_1, muon_sf_friends.MuonIsoSF_friends_1,],
    )

    configuration.add_outputs(["mt"], [q.id_wgt_mu_friend_1, q.iso_wgt_mu_friend_1,])

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()

