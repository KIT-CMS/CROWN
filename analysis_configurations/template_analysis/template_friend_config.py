from __future__ import annotations  # needed for type annotations in > python 3.7

from typing import List, Union

from .producers import pairselection as pairselection
from .producers import genparticles as genparticles
from .quantities import output as q
from code_generation.friend_trees import FriendTreeConfiguration
from code_generation.rules import RemoveProducer


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
        ["mm"],
        [
            pairselection.LVMu1_friend,
            pairselection.LVMu2_friend,
            genparticles.LVGenParticle1_friend,
            genparticles.LVGenParticle2_friend,
        ],
    )

    configuration.add_outputs(
        ["mm"],
        [
            q.p4_1,
            q.p4_2,
            q.gen_p4_1,
            q.gen_p4_2,
        ],
    )

    configuration.add_modification_rule(
        "mm",
        RemoveProducer(
            producers=[
                genparticles.LVGenParticle1_friend,
                genparticles.LVGenParticle2_friend,
            ],
            samples=["data"],
        ),
    )

    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()
