from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
from typing import Dict, List, Set, Union

log = logging.getLogger(__name__)


class Quantity:
    def __init__(self, name: str):
        self.name = name
        self.shifts: Dict[str, Set[str]] = {}
        self.ignored_shifts: Dict[str, Set[str]] = {}
        self.children: Dict[str, List[Quantity]] = {}
        self.defined_for_scopes: List[str] = []
        log.debug("Setting up new Quantity {}".format(self.name))

    def __str__(self) -> str:
        return self.name

    def __repr__(self) -> str:
        return self.name

    def __lt__(self, other: Quantity) -> bool:
        return self.name < other.name

    def add(self, name: str) -> None:
        """
        Function is not used for a base quantity
        """
        pass

    def reserve_scope(self, scope: str) -> None:
        """
        Function to reserve a scope for a given quantity. The scopes, in which a quantity is used
        as an output are tracked in the output_scopes list.
        If a quantity is already used within a given scope as output, this will result in an exception.
        This check is triggered for every Producer.
        """
        log.debug("Checking {} / scope {}".format(self.name, scope))
        if (scope == "global" and self.defined_for_scopes == []) or (
            scope != "global" and scope not in self.defined_for_scopes
        ):
            self.defined_for_scopes.append(scope)
        else:
            log.error(
                "Quantity {} is already defined in {} scope !".format(self.name, scope)
            )
            raise Exception

    def get_leaf(self, shift: str, scope: str) -> str:
        """
        Function to get the leaf of a given shift within a given scope.
        A leaf is the name of the quantity used for that scope/shift combination.
        If no shift is defined, the name of the quantity is returned.

        Args:
            shift (str): Name of the shift for which the leaf should be returned
            scope (str): Scope for which the leaf should be returned
        Returns:
            str. Name of the leaf
        """
        if shift in self.get_shifts(scope):
            return self.name + shift
        return self.name

    def get_leaves_of_scope(self, scope: str) -> List[str]:
        """
        Function returns a list of all leaves, which are defined for a given scope.

        Args:
            scope (str): Scope for which the leaves should be returned
        Returns:
            list. List of leaves
        """
        result = [self.name] + [
            self.get_leaf(shift, scope) for shift in self.get_shifts(scope)
        ]
        return result

    def shift(self, name: str, scope: str) -> None:
        """
        Function to define a shift for a given scope. If the shift is marked as ignored, nothing will be added.
        When a new shift is defined, all child quantities of a given quantity will be shifted as well.

        Args:
            name (str): Name of the shift
            scope (str): Scope for which the shift should be defined
        Returns:
            None
        """
        if scope in self.ignored_shifts.keys():
            if name in self.ignored_shifts[scope]:
                log.debug("Ignoring shift {} for quantity {}".format(name, self.name))
                return
        log.debug("Adding shift {} to quantity {}".format(name, self.name))
        if scope not in self.shifts.keys():
            self.shifts[scope] = set()
        if name not in self.shifts[scope]:
            self.shifts[scope].add(name)
            if scope == "global":  # shift children in all scopes if scope is global
                for any_scope in self.children:
                    for c in self.children[any_scope]:
                        c.shift(name, any_scope)
            else:
                if scope in self.children.keys():
                    for c in self.children[scope]:
                        c.shift(name, scope)

    def ignore_shift(self, name: str, scope: str) -> None:
        """
        Function to ignore a shift for a given scope.

        Args:
            name (str): Name of the shift to be ignored
            scope (str): Scope for which the shift should be ignored
        Returns:
            None
        """
        log.debug("Make quantity {} ignore shift {}".format(self.name, name))
        if scope not in self.ignored_shifts.keys():
            self.ignored_shifts[scope] = set()
        self.ignored_shifts[scope].add(name)

    def copy(self, name: str) -> Quantity:
        """
        Generate a copy of the current quantity with a new name.

        Args:
            name (str): Name of the new quantity
        Returns:
            Quantity. a new Quantity object.
        """
        copy = Quantity(name)
        copy.shifts = self.shifts
        copy.children = self.children
        copy.ignored_shifts = self.ignored_shifts
        copy.defined_for_scopes = self.defined_for_scopes
        return copy

    def adopt(self, child: Quantity, scope: str) -> None:
        """
        Adopt a child quantity to the current quantity.
        An adopted quantity will inherit all shifts of the partent quantity.

        Args:
            child (Quantity): The child quantity
            scope (str): Scope for which the child should be adopted
        Returns:
            None
        """
        log.debug(
            "Adopting child quantity {} to quantity {}".format(child.name, self.name)
        )
        if scope not in self.children.keys():
            self.children[scope] = []
        self.children[scope].append(child)

    def get_shifts(self, scope: str) -> List[str]:
        """
        Function returns a list of all shifts, which are defined for a given scope.

        Args:
            scope (str): Scope for which shifts should be returned
        Returns:
            list: List of all shifts, which are defined for a given scope.
        """
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
    """
    A Quantity Group is a group of quantities, that all have the same settings, but different names.
    """

    def __init__(self, name: str):
        super().__init__(name)
        self.quantities: List[Quantity] = []
        self.vec_config: str = ""

    def set_vec_config(self, vec_config: str) -> None:
        """
        Function to set the vec config key

        Args:
            vec_config (str): Name of the vec config key
        Returns:
            None
        """
        self.vec_config = vec_config

    def copy(self, name: str) -> Quantity:
        """
        Copy is not allowed for Quantity Groups.
        """
        log.error("Copy is not allowed for a Quantity Group !")
        raise Exception

    def add(self, name: str) -> None:
        """
        Function to add a new Quantity to the group. This quantity contains the identical shifts as the group itself

        Args:
            name (str): Name of the new Quantity
        Returns:
            None
        """
        if name not in [q.name for q in self.quantities]:
            quantity = Quantity(name)
            quantity.shifts = self.shifts
            quantity.children = self.children
            quantity.ignored_shifts = self.ignored_shifts
            quantity.defined_for_scopes = self.defined_for_scopes
            self.quantities.append(quantity)

    def get_leaves_of_scope(self, scope: str) -> List[str]:
        """
        Function returns a list of all leaves, which are defined for a given scope.
        This is an overload of the function used for the quantity class.

        For the writeout,  loop over all quantities in the group and
        return them all (plus their shifts) in a list.

        Args:
            scope (str): Scope for which the leaves should be returned
        Returns:
            list. List of leaves
        """
        output: List[str] = []
        for quantity in self.quantities:
            output.extend(quantity.get_leaves_of_scope(scope))
        return output


class NanoAODQuantity(Quantity):
    """
    A NanoAODQuantity is a quantity that comes directly from the NanoAOD file.
    Normally, these quantities are not suited to be used in the output ntuple,
    are therefore shielded from using them directly as a output.
    """

    def __init__(self, name: str):
        super().__init__(name)
        self.shifted_naming: Dict[str, str] = {}

    def reserve_scope(self, scope: str) -> None:
        """
        Function used to ensure, that NanoAOD quantities are not used as output
        in a producer. If this is attempted, an Exception will be raised.

        Args:
            scope (str): Scope for which the quantity should be reserved (not used)
        Returns:
            None
        """
        log.error(
            "Quantity {} is a NanoAOD quantity and cant be used as output !".format(
                self.name
            )
        )
        raise Exception

    def register_external_shift(
        self, shift_name: str, external_name: Union[str, NanoAODQuantity]
    ) -> None:
        """
        Function used to register a NanoAOD quantity as a shift of another quantity.
        Iif the shifted version of a quantity already exists in the input,
        this function can be used to register an
        branch from the input as a shifted version of a quantity

        Args:
            shift_name (str): Name of the shift
            external_name (str): Name of the shifted quantity in the NanoAOD
        Returns:
            None
        """
        if shift_name not in self.shifted_naming.keys():
            if isinstance(external_name, NanoAODQuantity):
                self.shifted_naming[shift_name] = str(external_name)
            else:
                self.shifted_naming[shift_name] = external_name
        for any_scope in self.children:
            for c in self.children[any_scope]:
                c.shift(shift_name, any_scope)

    def get_leaf(self, shift: str, scope: str) -> str:
        """
        Overloaded version for the `get_leaf` function to return the
        shifted versions of the quantity if needed.

        Args:
            shift (str): Name of the shift for which the leaf should be returned
            scope (str): Scope for which the leaf should be returned (not used)
        Returns:
            str. Name of the leaf
        """
        if shift in self.shifted_naming.keys():
            return self.shifted_naming[shift]
        return self.name

    def get_shifts(self, scope: str) -> List[str]:
        """
        Overloaded version of the `get_shifts` function

        Args:
            scope (str): Scope for which shifts should be returned (not used)
        Returns:
            list: List of all shifts
        """
        return list(self.shifted_naming.keys())


# Definitions for type annotations
QuantitiesInput = Union[
    Quantity, NanoAODQuantity, List[Union[Quantity, NanoAODQuantity]]
]
QuantitiesStore = Dict[str, Set[Quantity]]
