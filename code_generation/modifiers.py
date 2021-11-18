from __future__ import annotations  # needed for type annotations in > python 3.7
from typing import List, Union, Dict
from code_generation.exceptions import (
    SampleConfigurationError,
    EraConfigurationError,
)

ConfigurationParameters = Union[str, int, float, bool]

ModifierDict = Dict[
    str,
    Union[
        List[ConfigurationParameters],
        List[Dict[str, ConfigurationParameters]],
        Dict[str, ConfigurationParameters],
        ConfigurationParameters,
    ],
]
ModifierResolved = Union[
    List[ConfigurationParameters],
    Dict[str, ConfigurationParameters],
    List[Dict[str, ConfigurationParameters]],
    ConfigurationParameters,
    None,
]


class Modifier(object):
    def __init__(self, modifier_dict: ModifierDict):
        self.modifier_dict = modifier_dict

    def __str__(self) -> str:
        return "Modifier: {}".format(self.modifier_dict)

    def __repr__(self) -> str:
        return "Modifier: {}".format(self.modifier_dict)


class SampleModifier(Modifier):
    def __init__(
        self,
        modifier_dict: ModifierDict,
        default: Union[str, int, float, bool, None] = None,
    ):
        super(SampleModifier, self).__init__(modifier_dict)
        self.modifier_dict = modifier_dict
        self.default = default
        self.samples: List[str] = list(self.modifier_dict.keys())

    def apply(self, configstr: str) -> ModifierResolved:
        if configstr in self.samples:
            return self.modifier_dict[configstr]
        elif self.default is not None:
            return self.default
        else:
            raise SampleConfigurationError(configstr, self.samples)


class EraModifier(Modifier):
    def __init__(
        self,
        modifier_dict: ModifierDict,
        default: Union[str, int, float, bool, None] = None,
    ):
        super(EraModifier, self).__init__(modifier_dict)
        self.modifier_dict = modifier_dict
        self.default = default
        self.eras: List[str] = list(self.modifier_dict.keys())

    def apply(self, configstr: str) -> ModifierResolved:
        if configstr in self.eras:
            return self.modifier_dict[configstr]
        elif self.default is not None:
            return self.default
        else:
            raise EraConfigurationError(configstr, self.eras)
