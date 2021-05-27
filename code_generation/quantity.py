import logging

log = logging.getLogger(__name__)


class Quantity:
    def __init__(self, name):
        self.name = name
        self.shifts = {}
        self.children = {}
        self.output_scopes = []
        log.debug("Setting up new Quantity {}".format(self.name))

    # the scopes, in which a quantity is used as an output is tracked in the output_scopes list.
    # This check is triggered for every producer.
    # If a quantity is already used within a given scope as output, this will result in an exception.
    def check_scope(self, scope):
        log.debug("Checking {} / scope {}".format(self.name, scope))
        if scope not in self.output_scopes:
            self.output_scopes.append(scope)
        else:
            log.error(
                "Quantity {} is already used as output in {} scope !".format(
                    self.name, scope
                )
            )
            raise Exception

    def get_leaf(self, shift, scope):
        if shift in self.get_shifts(scope):
            return self.name + shift
        return self.name

    def get_leafs_of_scope(self, scope):
        return [self.name] + [self.name + shift for shift in self.get_shifts(scope)]

    def shift(self, name, scope):
        log.debug("Adding shift {} to quantity {}".format(name, self.name))
        if not scope in self.shifts.keys():
            self.shifts[scope] = set()
        if not name in self.shifts[scope]:
            self.shifts[scope].add(name)
            if scope == "global":  # shift children in all scopes if scope is global
                for any_scope in self.children.values():
                    for c in any_scope:
                        c.shift(name, scope)
            else:
                if scope in self.children.keys():
                    for c in self.children[scope]:
                        c.shift(name, scope)

    def copy(self, name):
        copy = Quantity(name)
        copy.shifts = self.shifts
        copy.children = self.children
        return copy

    def adopt(self, child, scope):
        log.debug(
            "Adopting child quantity {} to quantity {}".format(child.name, self.name)
        )
        if not scope in self.children.keys():
            self.children[scope] = []
        self.children[scope].append(child)

    def get_shifts(self, scope):
        if scope in self.shifts.keys():
            if "global" in self.shifts.keys():
                return list(self.shifts[scope].union(self.shifts["global"]))
            else:
                return list(self.shifts[scope])
        else:
            if "global" in self.shifts.keys():
                return list(self.shifts["global"])
            else:
                return []


class NanoAODQuantity(Quantity):
    def __init__(self, name):
        self.name = name
        super().__init__(name)

    # Quantities from the NanoAOD are not designed to be directly usable as output
    def check_scope(self, scope):
        log.error(
            "Quantity {} is a NanoAOD quantity and cant be used as output !".format(
                self.name
            )
        )
        raise Exception
