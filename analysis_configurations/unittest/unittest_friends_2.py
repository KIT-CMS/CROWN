from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List, Union
from .producers import muon_sf_friends as muon_sf_friends
from .producers import pairquantities as pairquantities
from .quantities import output as q

# from code_generation.configuration import Cnofiguration
from code_generation.friend_trees import FriendTreeConfiguration


def is_friend():
    return True


def build_config(
    era: str,
    sample: str,
    scopes: List[str],
    shifts: List[str],
    available_sample_types: List[str],
    available_eras: List[str],
    available_scopes: List[str],
    quantities_map: Union[str, None] = None,
):
    # for the test, we provide a quantities map
    configuration = FriendTreeConfiguration(
        era,
        sample,
        scopes,
        shifts,
        available_sample_types,
        available_eras,
        available_scopes,
        quantities_map,
    )

    configuration.add_producers(
        ["mt", "et", "tt", "em"],
        [
            pairquantities.FastMTTQuantities,
            muon_sf_friends.Rename_IDSF,
        ],
    )

    configuration.add_outputs(
        ["mt", "et", "tt", "em"],
        [
            q.m_fastmtt,
            q.pt_fastmtt,
            q.eta_fastmtt,
            q.phi_fastmtt,
            q.id_wgt_mu_friend_1_renamed,
        ],
    )

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()
