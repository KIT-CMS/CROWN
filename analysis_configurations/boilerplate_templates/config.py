from __future__ import annotations  # needed for type annotations in > python 3.7
from typing import List, Union
from code_generation.configuration import Configuration
from .producers import all_producers as p
from .quantities import output as q
from .quantities import nanoAOD as nanoAOD

def build_config(
    era: str,
    sample: str,
    scopes: List[str],
    shifts: List[str],
    available_sample_types: List[str],
    available_eras: List[str],
    available_scopes: List[str],
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



    #########################
    # Finalize and validate the configuration
    #########################
    configuration.optimize()
    configuration.validate()
    configuration.report()
    return configuration.expanded_configuration()