from __future__ import annotations  # needed for type annotations in > python 3.7
from typing import List, Union, Dict
from code_generation.exceptions import (
    SampleConfigurationError,
    EraConfigurationError,
)


class Modifier(object):
    def __init__(self, modifier_dict: Dict[str, Union[str, int, float, bool]]):
        self.modifier_dict = modifier_dict

    def apply(self, configstr: str) -> Dict[str, Union[str, int, float, bool]]:
        pass

    def __str__(self) -> str:
        return "Modifier: {}".format(self.modifier_dict)

    def __repr__(self) -> str:
        return "Modifier: {}".format(self.modifier_dict)


class SampleModifier(Modifier):
    def __init__(
        self,
        modifier_dict: Dict[str, Union[str, int, float, bool]],
        default: Union[bool, None] = None,
    ):
        super(SampleModifier, self).__init__(modifier_dict)
        self.modifier_dict = modifier_dict
        self.default = default
        self.samples: List[str] = list(self.modifier_dict.keys())

    def apply(self, configstr: str) -> Dict[str, Union[str, int, float, bool]]:
        if configstr in self.samples:
            return self.modifier_dict[configstr]
        elif self.default is not None:
            return self.default
        else:
            raise SampleConfigurationError(configstr, self.samples)


class EraModifier(Modifier):
    def __init__(
        self,
        modifier_dict: Dict[str, Union[str, int, float, bool]],
        default: Union[bool, None] = None,
    ):
        super(EraModifier, self).__init__(modifier_dict)
        self.modifier_dict = modifier_dict
        self.default = default
        self.eras = list(self.modifier_dict.keys())

    def apply(self, configstr: str) -> Dict[str, Union[str, int, float, bool]]:
        if configstr in self.eras:
            return self.modifier_dict[configstr]
        elif self.default is not None:
            return self.default
        else:
            raise EraConfigurationError(configstr, self.eras)
