from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
import ROOT
import json
import os
from time import time
from code_generation.configuration import Configuration
from typing import List, Union, Dict, Set

from code_generation.exceptions import (
    ConfigurationError,
    InvalidOutputError,
    InsufficientShiftInformationError,
)
from code_generation.producer import Producer, ProducerGroup
from code_generation.rules import ProducerRule

log = logging.getLogger(__name__)


class FriendTreeConfiguration(Configuration):
    """
    Configuration class for a FriendTree production with the CROWN framework.
    Based on the main Configuration class, but with a few modifications nessessary
    for a FriendTree configuration. The biggest differences are

    * the nominal version of quantities is optional and should only run if the user specifies it
    * no global scope is required
    * only one scope is allowed
    * The ordering is not optimized, but taken directly from the configuration file
    * information about the input file is required. This information can be provided by a json file, or by providing an input root file.
    * When using an input root file, only a single sample type and a single scope are allowed
    """

    def __init__(
        self,
        era: str,
        sample: str,
        scope: Union[str, List[str]],
        shifts: Set[str],
        available_sample_types: Union[str, List[str]],
        available_eras: Union[str, List[str]],
        available_scopes: Union[str, List[str]],
        input_information: Union[str, Dict[str, List[str]]],
    ):
        """Generate a configuration for a FriendTree production.

        Args:
            era (str): The era of the sample
            sample (str): The sample type
            scope (Union[str, List[str]]): The scope of the sample
            shifts (Set[str]): The shifts to be applied, can be "all", "none" or a list of shifts
            available_sample_types (Union[str, List[str]]): A list of available sample types
            available_eras (Union[str, List[str]]): A list of available eras
            available_scopes (Union[str, List[str]]): A list of available scopes
            input_information (Union[str, Dict[str, List[str]]]): Information about the input file. Can be a json file or a root file

        """  #

        super().__init__(
            era,
            sample,
            scope,
            shifts,
            available_sample_types,
            available_eras,
            available_scopes,
        )
        self.run_nominal = False
        # in the main constructor, the global scope is added to the scopes list.
        # This is not needed for a friend tree configuration, so we remove it again here
        if self.global_scope in self.scopes:
            self.scopes.remove(self.global_scope)
        self.global_scope = None
        log.warn(f"Selected scopes: {self.selected_scopes}")
        # if more than one scope is specified, raise an error
        if len(self.selected_scopes) > 1:
            raise ConfigurationError(
                f"FriendTree configurations can only have one scope, but multiple {self.selected_scopes} were specified"
            )

        self.input_quantities_mapping = self._readout_input_information(
            input_information
        )
        # all requested shifts are stored in a seperate varaiable,
        # they have to be added to all producers later
        self.requested_shifts = self._determine_requested_shifts(shifts)

    def _determine_requested_shifts(self, shiftset: Set[str]) -> Dict[str, List[str]]:
        """Determine the requested shifts from the user input

        Args:
            shifts (Union[str, List[str]]): User input for the shifts

        Returns:
            List[str]: List of requested shifts
        """
        requested_shifts: Dict[str, List[str]] = {}
        # first convert shifts to a list
        shifts = list(shiftset)
        testshifts = [shift.lower() for shift in shifts]
        # check if the user has specified "all" or "none"
        if "all" in testshifts or "none" in testshifts:
            if len(testshifts) > 1:
                raise ConfigurationError(
                    "When using 'all' or 'none' as a shift, no other shifts can be specified"
                )
            if testshifts[0] == "all":
                for scope in self.selected_scopes:
                    requested_shifts[scope] = list(
                        self.input_quantities_mapping[scope].keys()
                    )
                    # remove "" from the list
                    requested_shifts[scope].remove("")
                    requested_shifts[scope].append("nominal")
                    self.run_nominal = True
            elif testshifts[0] == "none":
                for scope in self.selected_scopes:
                    requested_shifts[scope] = ["nominal"]
                    self.run_nominal = True
        else:
            # in this case, the user has specified a list of shifts
            # we have to check if the shifts are valid
            for scope in self.selected_scopes:
                requested_shifts[scope] = []
                for shift in shifts:
                    if shift == "nominal":
                        self.run_nominal = True
                    elif shift not in self.input_quantities_mapping[scope].keys():
                        raise InsufficientShiftInformationError(
                            shift, list(self.input_quantities_mapping[scope].keys())
                        )
                    requested_shifts[scope].append(shift)
        return requested_shifts

    def _readout_input_information(
        self,
        input_information: Union[str, Dict[str, List[str]]],
    ) -> Dict[str, Dict[str, List[str]]]:
        # first check if the input is a root file or a json file
        data = {}
        if isinstance(input_information, str):
            if input_information.endswith(".root"):
                data = self._readout_input_root_file(input_information)
            elif input_information.endswith(".json"):
                data = self._readout_input_json_file(input_information)
            else:
                error_message = f"\n      Input information file {input_information} is not a json or root file \n"
                error_message += (
                    "      Did you forget to specify the input information file? \n"
                )
                error_message += "      The input information has to be a json file or a root file \n"
                error_message += "      and added to the cmake command via the -DQUANTITIESMAP=... option"
                raise ConfigurationError(error_message)
        return data

    def _readout_input_root_file(
        self, input_file: str
    ) -> Dict[str, Dict[str, List[str]]]:
        """Read the shift_quantities_map from the input root file and return it as a dictionary

        Args:
            input_file (str): Path to the input root file

        Returns:
            Dict[str, List[str]]: Dictionary containing the shift_quantities_map
        """

        data = {}
        if len(self.selected_scopes) > 1:
            raise ConfigurationError(
                "When using an input root file, only a single scope is possible"
            )

        start = time()
        log.debug(f"Reading quantities information from {input_file}")
        ROOT.gSystem.Load(os.path.abspath(__file__), "/maplib.so")
        f = ROOT.TFile.Open(input_file)
        name = "shift_quantities_map"
        m = f.Get(name)
        for shift, quantities in m:
            data[str(shift)] = [str(quantity) for quantity in quantities]
        f.Close()
        log.debug(
            f"Reading quantities information took {round(time() - start,2)} seconds"
        )
        return {list(self.selected_scopes)[0]: data}

    def _readout_input_json_file(
        self, input_file: str
    ) -> Dict[str, Dict[str, List[str]]]:
        """Read the shift_quantities_map from the input json file and return it as a dictionary

        Args:
            input_file (str): Path to the input json file

        Returns:
            Dict[str, List[str]]: Dictionary containing the shift_quantities_map
        """
        with open(input_file) as f:
            data = json.load(f)
        # json file structure is: {era: {sampletype: {scope: {shift: [quantities]}}}
        if self.era not in data:
            errorstring = (
                f"Era {self.era} not found in input information file {input_file}.\n"
            )
            errorstring += f"Available eras are: {data.keys()}"
            raise ConfigurationError(errorstring)
        if self.sample not in data[self.era].keys():
            errorstring = f"Sampletype {self.sample} not found in input information file {input_file}.\n"
            errorstring += f"Available sampletypes are: {data[self.era].keys()}"
            raise ConfigurationError(errorstring)
        if not set(self.selected_scopes).issubset(
            set(data[self.era][self.sample].keys())
        ):
            errorstring = f"Scopes {self.selected_scopes} not found in input information file {input_file}.\n"
            errorstring += f"Available scopes are: {data[self.era][self.sample].keys()}"
            raise ConfigurationError(errorstring)
        else:
            data = data[self.era][self.sample]
        return data

    def optimize(self) -> None:
        """
        Function used to optimize the FriendTreeConfiguration. In this case, no ordering
        optimization is performed. Optimizaion steps are:


            1. Apply rules
            2. Add all requested shifts to all producers. This addition is trivial, since
                the shifted quantities are already available in the input file
            3. Remove empty scopes

        Args:
            None

        Returns:
            None
        """
        self._apply_rules()
        self._add_requested_shifts()
        self._remove_empty_scopes()

    def _add_requested_shifts(self) -> None:
        # first shift the output quantities
        for scope in self.selected_scopes:
            for shift in self.requested_shifts[scope]:
                if shift != "nominal":
                    shiftname = "__" + shift
                    for producer in self.producers[scope]:
                        log.debug("Adding shift %s to producer %s", shift, producer)
                        producer.shift(shiftname, scope)
                        # second step is to shift the inputs of the producer
                        self._shift_producer_inputs(producer, shift, scope)
                        self.shifts[scope][shiftname] = {}

    def _shift_producer_inputs(self, producer, shift, scope):
        # if the producer is not of Type ProducerGroup we can directly shift the inputs
        if isinstance(producer, Producer):
            inputs = producer.get_inputs(scope)
            # only shift if necessary
            if shift in self.input_quantities_mapping[scope].keys():
                inputs_to_shift = []
                for input in inputs:
                    if input.name in self.input_quantities_mapping[scope][shift]:
                        inputs_to_shift.append(input)
                log.debug(f"Shifting inputs {inputs_to_shift} of producer {producer}")
                producer.shift_inputs("__" + shift, scope, inputs_to_shift)
        elif isinstance(producer, ProducerGroup):
            for producer in producer.producers:
                self._shift_producer_inputs(producer, shift, scope)

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
