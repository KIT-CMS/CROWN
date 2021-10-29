class ConfigurationError(Exception):
    """
    Exception raised when the configuration provided by the user is not valid.
    """

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class InvalidOutputError(ConfigurationError):
    """
    Exception raised when the list of output provided by the user is not valid.
    """

    def __init__(self, scope, outputs):
        self.message = "The required outputs {} for the scope '{}' are not provided by any producer \n Please modify the list of producers or the list of required outputs".format(
            outputs, scope
        )
        super().__init__(self.message)


class ChannelConfigurationError(ConfigurationError):
    """
    Exception raised when the channel configuration provided by the user is not valid.
    """

    def __init__(self, channels, available_channels):
        self.message = "Channel {} cannot be used in the configuration since it is not setup properly. Available channels are {}".format(
            channels, available_channels
        )
        super().__init__(self.message)


class SampleConfigurationError(ConfigurationError):
    """
    Exception raised when the sample type configuration provided by the user is not valid.
    """

    def __init__(self, samples, available_samples):
        self.message = "Sampletype {} cannot be used in the configuration since it is not setup properly. Available samples types are {}".format(
            samples, available_samples
        )
        super().__init__(self.message)


class EraConfigurationError(ConfigurationError):
    """
    Exception raised when the era configuration provided by the user is not valid.
    """

    def __init__(self, era, available_eras):
        self.message = "Era {} cannot be used in the configuration since it is not setup properly. Available eras are {}".format(
            era, available_eras
        )
        super().__init__(self.message)
