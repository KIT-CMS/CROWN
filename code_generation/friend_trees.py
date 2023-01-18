from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
from code_generation.configuration import Configuration
from typing import List, Union

from code_generation.exceptions import (
    ConfigurationError,
    InvalidOutputError,
    InsufficientShiftInformationError,
)
from code_generation.rules import ProducerRule

log = logging.getLogger(__name__)


class FriendTreeConfiguration(Configuration):
    """
    Configuration class for a FriendTree production with the CROWN framework.
    Based on the main Configuration class, but with a few modifications nessessary
    for a FriendTree configuration. The biggest differences
        * the nominal version of quantities is optional and should only run if the user specifies it
        * no global scope is required
        * The ordering is not optimized, but taken directly from the configuration file
    """

    def __init__(
        self,
        era: str,
        sample: str,
        scopes: Union[str, List[str]],
        shifts: Union[str, List[str]],
        available_sample_types: Union[str, List[str]],
        available_eras: Union[str, List[str]],
        available_scopes: Union[str, List[str]],
        run_nominal=False,
    ):
        super().__init__(
            era,
            sample,
            scopes,
            shifts,
            available_sample_types,
            available_eras,
            available_scopes,
        )
        self.run_nominal = run_nominal
        # in the main constructor, the global scope is added to the scopes list.
        # This is not needed for a friend tree configuration, so we remove it again here
        if self.global_scope in self.scopes:
            self.scopes.remove(self.global_scope)
        self.global_scope = None

        # catch the case where the user specifies All as shifts
        if shifts == "all":
            raise InsufficientShiftInformationError(shifts)

    def optimize(self) -> None:
        """
        Function used to optimize the FriendTreeConfiguration. In this case, no ordering optimization is performed. Optimizaion steps are:

            1. Remove empty scopes
            2. Apply rules

        Args:
            None

        Returns:
            None
        """
        self._apply_rules()
        self._remove_empty_scopes()

    def _validate_outputs(self) -> None:
        """
        Function used to validate the defined outputs. If an output is requested in the configuration,
        but is not available, since no producer will be able to produce it, an error is raised.

        Args:
            None

        Returns:
            None
        """
        for scope in [scope for scope in self.scopes]:
            required_outputs = set(output for output in self.outputs[scope])
            # merge the two sets of outputs
            provided_outputs = self.available_outputs[scope]
            missing_outputs = required_outputs - provided_outputs
            if len(missing_outputs) > 0:
                raise InvalidOutputError(scope, missing_outputs)

    def add_modification_rule(
        self, scopes: Union[str, List[str]], rule: ProducerRule
    ) -> None:
        """
        Function used to add a rule to the configuration.

        Args:
            scopes: The scopes to which the rule should be added. This can be a list of scopes or a single scope.
            rule: The rule to be added. This must be a ProducerRule object.

        Returns:
            None

        """
        if not isinstance(rule, ProducerRule):
            raise TypeError("Rule must be of type ProducerRule")
        if not isinstance(scopes, list):
            scopes = [scopes]
        rule.set_scopes(scopes)
        # TODO Check if this works without a global scope
        if self.global_scope is not None:
            rule.set_global_scope(self.global_scope)
        self.rules.add(rule)

    def expanded_configuration(self) -> Configuration:
        """Function used to generate an expanded version of the configuration, where all shifts are applied.
        This expanded configuration is used by the code generator to generate the C++ code.

        Returns:
            Configuration: Expanded configuration
        """
        expanded_configuration = {}
        for scope in self.scopes:
            expanded_configuration[scope] = {}
            if self.run_nominal:
                log.debug("Adding nominal in scope {}".format(scope))
                expanded_configuration[scope]["nominal"] = self.config_parameters[scope]
            if len(self.shifts[scope]) > 0:
                for shift in self.shifts[scope]:
                    log.debug("Adding shift {} in scope {}".format(shift, scope))
                    log.debug("  {}".format(self.shifts[scope][shift]))
                    try:
                        expanded_configuration[scope][shift] = (
                            self.config_parameters[scope] | self.shifts[scope][shift]
                        )
                    except KeyError:
                        expanded_configuration[scope][shift] = {}
                        expanded_configuration[scope][shift] = (
                            self.config_parameters[scope] | self.shifts[scope][shift]
                        )
            # log.warn("Expanded configuration for scope {}".format(scope))
            # log.warn(expanded_configuration[scope])
        # check if any shift (including the nominal) is run, if not, exit with an error
        if not any(
            [len(expanded_configuration[scope]) > 0 for scope in expanded_configuration]
        ):
            error_msg = "Nothing to run, is the configuration valid? \n Provided Configuration: \n {}".format(
                self
            )
            raise ConfigurationError(error_msg)

        self.config_parameters = expanded_configuration
        return self
