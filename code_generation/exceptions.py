from __future__ import annotations  # needed for type annotations in > python 3.7
from typing import List, Set, Union
from code_generation.quantity import Quantity


class ConfigurationError(Exception):
    """
    Exception raised when the configuration provided by the user is not valid.
    """

    def __init__(self, message: str):
        self.message = message
        super().__init__(self.message)


class InvalidOutputError(ConfigurationError):
    """
    Exception raised when the list of output provided by the user is not valid.
    """

    def __init__(self, scope: str, outputs: Union[Set[Quantity], List[Quantity]]):
        self.message = "The required outputs {} for the scope '{}' are not provided by any producer \n Please modify the list of producers or the list of required outputs".format(
            outputs, scope
        )
        super().__init__(self.message)


class ChannelConfigurationError(ConfigurationError):
    """
    Exception raised when the channel configuration provided by the user is not valid.
    """

    def __init__(
        self, channels: Set[str], available_channels: Union[Set[str], List[str]]
    ):
        self.message = "Channels {} cannot be used in the configuration since it is not setup properly. Available channels are {}".format(
            channels, available_channels
        )
        super().__init__(self.message)


class SampleConfigurationError(ConfigurationError):
    """
    Exception raised when the sample type configuration provided by the user is not valid.
    """

    def __init__(self, sample: str, available_samples: Union[Set[str], List[str]]):
        self.message = "Sampletype {} cannot be used in the configuration since it is not setup properly. Available samples types are {}".format(
            sample, available_samples
        )
        super().__init__(self.message)


class EraConfigurationError(ConfigurationError):
    """
    Exception raised when the era configuration provided by the user is not valid.
    """

    def __init__(self, era: str, available_eras: Union[Set[str], List[str]]):
        self.message = "Era {} cannot be used in the configuration since it is not setup properly. Available eras are {}".format(
            era, available_eras
        )
        super().__init__(self.message)


class InvalidProducerConfigurationError(ConfigurationError):
    """
    Exception raised when the producer configuration provided by the user is not valid.
    """

    def __init__(self, producer: str):
        self.message = "Producer {} is not setup properly".format(
            producer,
        )
        super().__init__(self.message)
