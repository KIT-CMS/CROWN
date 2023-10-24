from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List, Union
import os
from .producers import muon_sf_friends as muon_sf_friends
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
    quantities_map = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "dyjets_shift_quantities_map.json",
    )
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

    configuration.add_config_parameters(
        ["mt", "mm"],
        {
            "muon_sf_file": "data/embedding/muon_2018UL.json.gz",
            "muon_id_sf": "ID_pt_eta_bins",
            "muon_iso_sf": "Iso_pt_eta_bins",
        },
    )

    configuration.add_producers(
        ["mt", "mm"],
        [
            muon_sf_friends.MuonIDSF_friends_1,
            muon_sf_friends.MuonIsoSF_friends_1,
        ],
    )

    configuration.add_outputs(
        ["mt", "mm"],
        [
            q.id_wgt_mu_friend_1,
            q.iso_wgt_mu_friend_1,
        ],
    )

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()
