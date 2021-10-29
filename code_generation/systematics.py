from code_generation.producer import Producer
import logging

log = logging.getLogger(__name__)


class SystematicShift(object):
    def __init__(self, name, shift_config, producers, scopes=None, ignore_producers={}):
        self.name = "__" + name
        self.producers = producers
        self.ignore_producers = ignore_producers
        self.shift_config = shift_config
        self.scopes = self.determine_scopes(scopes)
        self.validate()

    def __str__(self) -> str:
        returnstr = "SystematicShift: {}\n".format(self.name)
        for scope in self.scopes:
            returnstr += "  Scope: {}\n".format(scope)
            returnstr += "  Producers: {}\n".format(self.producers[scope])
            returnstr += "  Ignore producers: {}\n".format(self.ignore_producers[scope])
            returnstr += "  Shift config: {}\n".format(self.shift_config[scope])
            returnstr += "  ------------ \n"
        return returnstr

    def __repr__(self) -> str:
        returnstr = "SystematicShift: {}\n".format(self.name)
        for scope in self.scopes:
            returnstr += "  Scope: {}\n".format(scope)
            returnstr += "  Producers: {}\n".format(self.producers[scope])
            returnstr += "  Ignore producers: {}\n".format(self.ignore_producers[scope])
            returnstr += "  Shift config: {}\n".format(self.shift_config[scope])
            returnstr += "  ------------ \n"
        return returnstr

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

    def convert_to_dict(self, producers, scopes):
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
                    log.warning(
                        "Setting empy Producer set for shift {} in scope {}".format(
                            self.name, scope
                        )
                    )
        if not isinstance(producers, dict):
            raise ValueError(
                "SystematicShift producers must be a dict, list or a single Producer"
            )
        return producer_dict

    def validate(self):
        # if isinstance(self.scopes, str):
        #     self.scopes = [self.scopes]
        # if not isinstance(self.scopes, list):
        #     raise ValueError("SystematicShift scopes must be a list or single string")
        self.producers = self.convert_to_dict(self.producers, self.scopes)
        self.ignore_producers = self.convert_to_dict(self.ignore_producers, self.scopes)
        unconfigured_scopes = set(self.scopes) - set(self.shift_config.keys())
        # for uncongured scopes, set an empty config
        if len(unconfigured_scopes) > 0:
            for scope in unconfigured_scopes:
                self.shift_config[scope] = {}
            # raise ValueError(
            #     "Shiftconfig for scopes {} not found in shift dict".format(
            #         unconfigured_scopes
            #     )
            # )

    def name(self):
        return self.name

    def add_producer(self, producer, scopes=None):
        if scopes is None:
            scopes = self.scopes
        if isinstance(scopes, str):
            scopes = [scopes]
        for scope in scopes:
            self.producers[scope].append(producer)
        self.validate()

    def add_ignore_producer(self, producer, scopes=None):
        if scopes is None:
            scopes = self.scopes
        if isinstance(scopes, str):
            scopes = [scopes]
        for scope in scopes:
            self.ignore_producers[scope].append(producer)
        self.validate()

    def add_scope(self, scope):
        self.scopes.append(scope)
        self.producers[scope] = []
        self.ignore_producers[scope] = []
        self.validate()

    def add_config(self, config, scope):
        self.shift_config[scope] = config

    def get_scopes(self):
        return self.scopes

    def get_shift_config(self, scope):
        return self.shift_config[scope]

    def apply(self, scope):
        log.warning("Applying systematic shift \n{}".format(self))
        if scope in self.scopes:
            for producer in self.ignore_producers[scope]:
                producer.ignore_shift(self.name, scope)
            for producer in self.producers[scope]:
                producer.shift(self.name, scope)


# # Function for introducing systematic variations to producers and depending quantities
# def AddSystematicShift(
#     config, name, change_dict, base_producers, sanitize_producers=[]
# ):
#     name = "__" + name
#     shift_config = copy.deepcopy(config[""])
#     for scope in change_dict:
#         shift_config[scope].update(change_dict[scope])
#     config[name] = shift_config
#     for producer in sanitize_producers:
#         producer[0].ignore_shift(name, producer[1])
#     for producer in base_producers:
#         producer[0].shift(name, producer[1])


# # Function for introducing systematic variations to producers and depending quantities by adding an already shifted input quantity
# def SystematicShiftByInputQuantity(config, shiftname, external_dict):
#     shiftname = "__" + shiftname
#     config[shiftname] = copy.deepcopy(config[""])
#     for quantity in external_dict.keys():
#         quantity.register_external_shift(
#             shift_name=shiftname,
#             external_name=external_dict[quantity],
#         )
