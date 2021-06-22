import logging

log = logging.getLogger(__name__)


class Quantity:
    def __init__(self, name):
        self.name = name
        self.shifts = {}
        self.ignored_shifts = {}
        self.children = {}
        self.defined_for_scopes = []
        log.debug("Setting up new Quantity {}".format(self.name))

    # the scopes, in which a quantity is used as an output is tracked in the output_scopes list.
    # This check is triggered for every producer.
    # If a quantity is already used within a given scope as output, this will result in an exception.
    def reserve_scope(self, scope):
        log.debug("Checking {} / scope {}".format(self.name, scope))
        if scope not in self.defined_for_scopes:
            self.defined_for_scopes.append(scope)
        else:
            log.error(
                "Quantity {} is already defined in {} scope !".format(self.name, scope)
            )
            raise Exception

    def get_leaf(self, shift, scope):
        if shift in self.get_shifts(scope):
            return self.name + shift
        return self.name

    def get_leafs_of_scope(self, scope):
        return [self.name] + [self.name + shift for shift in self.get_shifts(scope)]

    def shift(self, name, scope):
        if scope in self.ignored_shifts.keys():
            if name in self.ignored_shifts[scope]:
                log.debug("Ignoring shift {} for quantity {}".format(name, self.name))
                return
        log.debug("Adding shift {} to quantity {}".format(name, self.name))
        if not scope in self.shifts.keys():
            self.shifts[scope] = set()
        if not name in self.shifts[scope]:
            self.shifts[scope].add(name)
            if scope == "global":  # shift children in all scopes if scope is global
                for any_scope in self.children:
                    for c in self.children[any_scope]:
                        c.shift(name, any_scope)
            else:
                if scope in self.children.keys():
                    for c in self.children[scope]:
                        c.shift(name, scope)

    def ignore_shift(self, name, scope):
        log.debug("Make quantity {} ignore shift {}".format(self.name, name))
        if not scope in self.ignored_shifts.keys():
            self.ignored_shifts[scope] = set()
        self.ignored_shifts[scope].add(name)

    def copy(self, name):
        copy = Quantity(name)
        copy.shifts = self.shifts
        copy.children = self.children
        copy.ignored_shifts = self.ignored_shifts
        copy.defined_for_scopes = self.defined_for_scopes
        return copy

    def adopt(self, child, scope):
        log.debug(
            "Adopting child quantity {} to quantity {}".format(child.name, self.name)
        )
        if not scope in self.children.keys():
            self.children[scope] = []
        self.children[scope].append(child)

    def get_shifts(self, scope):
        if "global" in self.shifts.keys():
            if scope != "global" and scope in self.shifts.keys():
                log.error(
                    "Quantity {} has shifts in global and {}. Something must be broken!".format(
                        self.name, scope
                    )
                )
                raise Exception
            return list(self.shifts["global"])
        elif scope in self.shifts.keys():
            return list(self.shifts[scope])
        else:
            return []


class QuantityGroup(Quantity):
    # A Quantity Group is a group of quantities, that all have the same settings, but different names.
    def __init__(self, name):
        super().__init__(name)
        self.quantities = []

    def copy(self, name):
        log.error("Copy is not allowed for a Quantity Group !")
        raise Exception

    def add(self, name):
        # add a new Quantity to the group. This quantity contains the identical shifts as the group itself
        quantity = Quantity(name)
        quantity.shifts = self.shifts
        quantity.children = self.children
        quantity.ignored_shifts = self.ignored_shifts
        quantity.defined_for_scopes = self.defined_for_scopes
        self.quantities.append(quantity)

    def get_leafs_of_scope(self, scope):
        # For the writeout, we have to loop over all quantities in the group and return them all (plus their shifts) in a list
        output = []
        for quantity in self.quantities:
            output.extend(
                [quantity.name]
                + [quantity.name + shift for shift in quantity.get_shifts(scope)]
            )
        return output


class NanoAODQuantity(Quantity):
    def __init__(self, name):
        super().__init__(name)

    # Quantities from the NanoAOD are not designed to be directly usable as output
    def reserve_scope(self, scope):
        log.error(
            "Quantity {} is a NanoAOD quantity and cant be used as output !".format(
                self.name
            )
        )
        raise Exception
