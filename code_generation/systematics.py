from code_generation.producer import Producer
import logging

log = logging.getLogger(__name__)

"""
    Class containing systematic shifts. A systematic shift is a variation of a 
    set of configuration parameters. A dummy shift looks like this:
    SystematicShift(
            name="shiftname",
            shift_config={
                ("scope1", "scope2"): {
                    "paramter1": True,
                    "paramter2": 42.0,
                },
            },
            producers={
                "scope1": modified_producer1,
                "scope2": modified_producer2,
            },
            ignore_producers={
                "scope1": ignored_producer1,
                "scope2": ignored_producer2,
            },
        )

    Args:
        name (str): Name of the systematic shift.

        shift_config (dict): Dictionary containing the configuration parameters. The dictionary keys have to be strings, or tuples, in case the same configuration change is used for multiple scopes. The dictionary values have to be dictionaries containing the changed configuration parameters.

        producers (dict): Dictionary containing the producers that are affected by the systematic shift. The dictionary keys have to be strings, or tuples, in case the same producer is affected by the systematic shift in multiple scopes. The dictionary values have to be the modified producers.

        ignore_producers (dict): Dictionary containing the producers that are not affected by the systematic shift. The dictionary keys have to be strings, or tuples, in case the same producer is not affected by the systematic shift in multiple scopes. The dictionary values have to be the ignored producers.

"""
class SystematicShift(object):
    def __init__(self, name, shift_config, producers, scopes=None, ignore_producers={}):
        self.shiftname = "__" + name
        self.producers = self.expand_dict(producers)
        self.ignore_producers = self.expand_dict(ignore_producers)
        self.shift_config = self.expand_dict(shift_config)
        self.scopes = self.determine_scopes(scopes)
        self.validate()
    
    """
    Function used to expand dictionaries. If the key is a string, it is returned as is. If the key is a tuple, the tuple is expanded and the resulting dict with one key per tuple entry is returned.

    Args:
        dict_to_expand (dict): Dictionary to expand.

    Returns:
        dict: Expanded dictionary.
    """

    def expand_dict(self, dict_to_expand) -> dict:
        exanded_dict = {}
        for key in dict_to_expand:
            if isinstance(key, str):
                exanded_dict[key] = dict_to_expand[key]
            elif isinstance(key, tuple):
                for x in key:
                    exanded_dict[x] = dict_to_expand[key]
            else:
                raise ValueError("Dict key must be a string or tuple")
        return exanded_dict

    def __str__(self) -> str:
        returnstr = "SystematicShift: {}\n".format(self.shiftname)
        for scope in self.scopes:
            returnstr += "  Scope: {}\n".format(scope)
            returnstr += "  Producers: {}\n".format(self.producers[scope])
            returnstr += "  Ignore producers: {}\n".format(self.ignore_producers[scope])
            returnstr += "  Shift config: {}\n".format(self.shift_config[scope])
            returnstr += "  ------------ \n"
        return returnstr

    def __repr__(self) -> str:
        returnstr = "SystematicShift: {}\n".format(self.shiftname)
        for scope in self.scopes:
            returnstr += "  Scope: {}\n".format(scope)
            returnstr += "  Producers: {}\n".format(self.producers[scope])
            returnstr += "  Ignore producers: {}\n".format(self.ignore_producers[scope])
            returnstr += "  Shift config: {}\n".format(self.shift_config[scope])
            returnstr += "  ------------ \n"
        return returnstr
    
    """
    Function used to determine the scopes that are affected by the systematic shift.
    If no scope is specified by the user, the scopes are determined from the shift_config, producers and ignore_producers.

    Args:
        scopes (set): Set of scopes that are affected by the systematic shift.
    
    Returns:
        set: Set of scopes that are affected by the systematic shift.
    """

    def determine_scopes(self, scopes) -> set:
        if scopes is None:
            scopes = (
                set(self.shift_config.keys())
                | set(self.producers.keys())
                | set(self.ignore_producers.keys())
            )
        if isinstance(scopes, str):
            scopes = [scopes]
        if isinstance(scopes, list):
            scopes = set(scopes)
        return scopes

    """
    
    Helper function used in the validation of the systematic shift. If makes sure, that the provided dict of producers is in a defined format, meaning a dict with one key per scope and a list of producers as value.

    Args:
        producers: Producers can be a single Producer object, a list of Producer objects or a dict with one key per scope and a list of producers as value.

    Returns:
        dict: Dictionary containing the producers with one key per scope and a list of producers as value.

    """

    def convert_to_dict(self, producers, scopes) -> dict:
        producer_dict = {}
        if isinstance(self.producers, Producer):
            for scope in scopes:
                producer_dict[scope] = [producers]
        elif isinstance(producers, list):
            for scope in scopes:
                producer_dict[scope] = [producer for producer in producers]
        elif isinstance(producers, dict):
            for scope in scopes:
                try:
                    if isinstance(producers[scope], Producer):
                        producer_dict[scope] = [producers[scope]]
                    elif isinstance(producers[scope], list):
                        producer_dict[scope] = producers[scope]
                    else:
                        raise ValueError(
                            "Producer must be a list or Producer for {}".format(self)
                        )
                except KeyError:
                    producer_dict[scope] = []
                    log.debug(
                        "Setting empy Producer set for shift {} in scope {}".format(
                            self.shiftname, scope
                        )
                    )
        if not isinstance(producers, dict):
            raise ValueError(
                "SystematicShift producers must be a dict, list or a single Producer"
            )
        return producer_dict
    
    """
    Function used to validate the systematic shift. The provided producers and ignore_producers are converted into a defined format. If no configuration is provided for a scope, an empty dict is used.

    Args:
        None

    Returns:
        None
    """

    def validate(self) -> None:
        self.producers = self.convert_to_dict(self.producers, self.scopes)
        self.ignore_producers = self.convert_to_dict(self.ignore_producers, self.scopes)
        unconfigured_scopes = set(self.scopes) - set(self.shift_config.keys())
        # for uncongured scopes, set an empty config
        if len(unconfigured_scopes) > 0:
            for scope in unconfigured_scopes:
                self.shift_config[scope] = {}

    """
    Function used to add a producer to the list of producers affected by the systematic shift. After adding the producer, the shift is validated.

    Args:
        producer: Producer to add.
        scope: Scope to which the producer should be added. If no sopce is provided, the producer is added to all scopes of the systematic.

    Returns:
        None
    """

    def add_producer(self, producer, scopes=None) -> None:
        if scopes is None:
            scopes = self.scopes
        if isinstance(scopes, str):
            scopes = [scopes]
        for scope in scopes:
            self.producers[scope].append(producer)
        self.validate()

     """
    Function used to add an ignored producer to the list of ignored producers, which are untouched by the systematic shift. After adding the ignored producer, the shift is validated.

    Args:
        producer: Producer to ignore.
        scope: Scope to which the ignored producer should be added. If no sopce is provided, the ignored producer is added to all scopes of the systematic.

    Returns:
        None
    """
    def add_ignore_producer(self, producer, scopes=None) -> None:
        if scopes is None:
            scopes = self.scopes
        if isinstance(scopes, str):
            scopes = [scopes]
        for scope in scopes:
            self.ignore_producers[scope].append(producer)
        self.validate()

    """
    Function used to add a scope to the list of scopes affected by the systematic shift. After adding the scope, the shift is validated.

    Args:
        scope: Scope to add.

    Returns:
        None
    """
    def add_scope(self, scope) -> None:
        self.scopes.add(scope)
        self.producers[scope] = []
        self.ignore_producers[scope] = []
        self.validate()

    """
    Function used to add a configuration to the list of configurations for a scope. After adding the configuration, the shift is validated.

    Args:
        config: Configuration to add.
        scope: Scope to which the configuration should be added.

    Returns:
        None
    """

    def add_config(self, config, scope) -> None:
        self.shift_config[scope] = config
        self.validate()

    """
    This function returns a set of all scopes that are affected by the systematic shift.

    Args:
        None

    Returns:
        set: Set of scopes that are affected by the systematic shift. 
    """

    def get_scopes(self) -> set:
        return self.scopes

    """
    This function returns the configuration for a given scope.

    Args:
        scope: Scope for which the configuration should be returned.

    Returns:
        dict: Configuration for the given scope.
    """  

    def get_shift_config(self, scope) -> dict:
        return self.shift_config[scope]


    """
    Function used to apply the systematic shift to the given producers. For the given scope, all producers aer shifted using producer.shift, while, for all ignored producers, the producer.ignore_shift function is called. If the scope is not defined in the shift, no shift is applied.

    Args:
        scope: Scope for which the shift should be applied.
    
    Returns:
        None
    """
    def apply(self, scope) -> None:
        log.debug("Applying systematic shift \n{}".format(self))
        if scope in self.scopes:
            for producer in self.ignore_producers[scope]:
                producer.ignore_shift(self.shiftname, scope)
            for producer in self.producers[scope]:
                producer.shift(self.shiftname, scope)
