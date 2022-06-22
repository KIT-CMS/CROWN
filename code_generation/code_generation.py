from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
from typing import Any, Dict, List, Union, Tuple
import os
import filecmp
from git import Repo, InvalidGitRepositoryError, NoSuchPathError

from code_generation.producer import SafeDict, Producer, ProducerGroup

from code_generation.configuration import Configuration

log = logging.getLogger(__name__)


class CodeSubset(object):

    """
    Class used to generate code for a smaller subset. For each subset, a new object must be created.

    Args:
        file_name: The name of the file to be generated.
        template: The template to be used for the generation.
        producer: The producer, of which the code will be generated.
        scope: The scope of the code generation.
        folder: The folder in which the code will be generated.
        parameters: The parameters to be used for the generation.

    Returns:
        None
    """

    def __init__(
        self,
        file_name: str,
        template: str,
        producer: Union[Producer, ProducerGroup],
        scope: str,
        folder: str,
        configuration_parameters: Dict[str, Any],
    ):
        self.file_name = file_name
        self.template = template
        self.producer = producer
        self.scope = scope
        self.name = producer.name + "_" + scope
        self.configuration_parameters = configuration_parameters
        self.count = 0
        self.folder = folder
        self.commands: List[str] = []
        self.headerfile = os.path.join(
            self.folder, "include", self.scope, "{}.hxx".format(self.file_name)
        )
        self.sourcefile = os.path.join(
            self.folder, "src", self.scope, "{}.cxx".format(self.file_name)
        )

    def create(self):
        """
        Create the code subset. Calls the writecalls function of the producer to generate the code.

        Args:
            None

        Returns:
            None
        """
        log.debug("Creating code subset {}".format(self.name))
        log.debug("Producer: {}".format(self.producer.name))
        log.debug("Scope: {}".format(self.scope))
        self.producer.reserve_output(self.scope)
        # create the function calls for the producer
        for call in self.producer.writecalls(self.configuration_parameters, self.scope):
            log.debug("Adding call for {}".format(self.name))
            log.debug("Call: {}".format(call))
            expanded_call = call.format_map(
                SafeDict(
                    {
                        "df": "df{}".format(self.count),
                        "vec_open": "{",
                        "vec_close": "}",
                    }
                )
            )
            self.commands.append(
                "    auto df{} = {};\n".format(self.count + 1, expanded_call)
            )
            self.count += 1
            log.debug("|---> {}".format(self.commands))
        self.commands.append("    return df{};\n".format(self.count))

    def write(self):
        """
        Write the code subset to a file, both the header and the source. Before writing the files,
        check if they already exists, and if they exist and are not different, skip writing them.
        This is to avoid unnecessary recompilation, since the compiler will check the timestamps of the files.

        Args:
            None

        Returns:
            None
        """
        log.debug("Writing code subset {}".format(self.name))
        log.debug("folder: {}, file_name: {}".format(self.folder, self.file_name))
        # write the header file if it does not exist or is different
        with open(self.headerfile + ".new", "w") as f:
            f.write(f"ROOT::RDF::RNode {self.name}(ROOT::RDF::RNode df);")
        if os.path.isfile(self.headerfile):
            if filecmp.cmp(self.headerfile + ".new", self.headerfile):
                log.debug("--> Identical header file, skipping")
                os.remove(self.headerfile + ".new")
            else:
                os.rename(self.headerfile + ".new", self.headerfile)
        else:
            os.rename(self.headerfile + ".new", self.headerfile)
        # write the source file if it does not exist or is different
        with open(self.sourcefile + ".new", "w") as f:
            commandstring = "".join(self.commands)
            f.write(
                self.template.replace("//    { commands }", commandstring).replace(
                    "{subsetname}", self.name
                )
            )
        if os.path.isfile(self.sourcefile):
            if filecmp.cmp(self.sourcefile + ".new", self.sourcefile):
                os.remove(self.sourcefile + ".new")
                log.debug("--> Identical source file, skipping")
            else:
                os.rename(self.sourcefile + ".new", self.sourcefile)
        else:
            os.rename(self.sourcefile + ".new", self.sourcefile)

    def call(self, inputscope: str, outputscope: str) -> str:
        """
        Return the call to the code subset. This call is used in the generated code of the executalbe.

        Args:
            inputscope: The scope of the input dataframe.
            outputscope: The scope of the output dataframe.

        Returns:
            str: the call to the code subset
        """
        call = f"    auto {outputscope} = {self.name}({inputscope}); \n"
        return call

    def include(self) -> str:
        """
        Return the include statement for the code subset

        Args:
            None

        Returns:
            str: the include statement for the code subset
        """
        return f'#include "{self.headerfile}"\n'


class CodeGenerator(object):
    """
    Class used to generate code from a given Configuration. The code is generated in a folder, which is the name of the executable.
    Inside the folder the source file for the executable is generated, as well as and include and source dir. Within those two folders,
    a subfolder for each scope is generated and within those, the code for each producer is generated. Each file contains all calls for one producer from the config.

    Args:
        main_template_path: the path to the cxx template for the executable
        sub_template_path: the path to the cxx template for the code subsets
        configuration: the configuration to generate code from
        analysis_name: the name of the analysis
        executable_name: the name of the executable
        output_folder: the folder to write the code to

    Returns:
        None
    """

    def __init__(
        self,
        main_template_path: str,
        sub_template_path: str,
        configuration: Configuration,
        analysis_name: str,
        executable_name: str,
        output_folder: str,
        threads: int = 1,
    ):
        self.main_template = self.load_template(main_template_path)
        self.subset_template = self.load_template(sub_template_path)
        self.configuration = configuration
        self.scopes = self.configuration.scopes
        self.outputs = self.configuration.outputs
        self.global_scope = self.configuration.global_scope
        self.executable_name = executable_name
        self.analysis_name = analysis_name
        self.output_folder = output_folder
        self.executable = os.path.join(
            output_folder,
            self.executable_name + "_generated_code",
            self.executable_name + ".cxx",
        )
        self.debug = False
        self._outputfiles_generated: Dict[str, str] = {}
        self.threads = threads
        self.subset_includes: List[str] = []
        self.output_commands: Dict[str, List[str]] = {}
        self.subset_calls: Dict[str, List[str]] = {}
        self.main_counter: Dict[str, int] = {}
        self.number_of_defines = 0
        self.number_of_outputs = 0
        for scope in self.scopes:
            self.main_counter[scope] = 0
            self.subset_calls[scope] = []
            self.output_commands[scope] = []
        # get git status of the repo
        try:
            repo = Repo("../../CROWN")
            self.commit_hash = repo.head.commit
            self.setup_is_clean = "false" if repo.is_dirty() else "true"
        except (ValueError, InvalidGitRepositoryError, NoSuchPathError):
            self.commit_hash = "undefined"
            self.setup_is_clean = "false"
        log.info("Code generator initialized")

    def generate_code(self) -> None:
        """
        Generate the code from the configuration and create the subsets. Run through the whole configuration and create a subset for each producer within the configuration.

        Start with the global scope and then all other scopes. All generated code is stored in the folder self.output_folder.

        Args:
            None

        Returns:
            None
        """
        # start with the global scope

        for subfolder in ["src", "include"]:
            for scope in self.scopes:
                if not os.path.exists(
                    os.path.join(
                        self.output_folder,
                        self.executable_name + "_generated_code",
                        subfolder,
                        scope,
                    )
                ):
                    os.makedirs(
                        os.path.join(
                            self.output_folder,
                            self.executable_name + "_generated_code",
                            subfolder,
                            scope,
                        )
                    )
        # self.generate_subsets(self.global_scope)
        for scope in self.scopes:
            self.generate_subsets(scope)

        calls, includes = self.generate_main_code()
        run_commands = self.generate_run_commands()

        self.write_code(calls, includes, run_commands)

    def load_template(self, template_path: str) -> str:
        """
        Load the template from the given path
        Args:
            template_path: the path to the template
        Returns:
            str - the template
        """
        with open(template_path, "r") as template_file:
            template = template_file.read()
        return template

    def write_code(self, calls: str, includes: str, run_commands: str) -> None:
        """
        Write the code of the main executable to the output folder

        Args:
            calls: the main calls
            includes: the includes
            run_commands: the run commands

        Returns:
            None
        """
        if self.threads > 1:
            log.info(f"Using {self.threads} threads for the executable")
            threadcall = "ROOT::EnableImplicitMT({});".format(self.threads)
        else:
            threadcall = ""
        with open(self.executable, "w") as f:
            f.write(
                self.main_template.replace("    // {CODE_GENERATION}", calls)
                .replace("// {INCLUDES}", includes)
                .replace("    // {RUN_COMMANDS}", run_commands)
                .replace("// {MULTITHREADING}", threadcall)
                .replace("// {DEBUGLEVEL}", self.set_debug_flag())
                .replace("{ERATAG}", '"Era={}"'.format(self.configuration.era))
                .replace(
                    "{SAMPLETAG}", '"Samplegroup={}"'.format(self.configuration.sample)
                )
                .replace("{ANALYSISTAG}", '"Analysis={}"'.format(self.analysis_name))
                .replace("{PROGRESS_CALLBACK}", self.set_process_tracking())
                .replace("{OUTPUT_QUANTITIES}", self.set_output_quantities())
                .replace("{SYSTEMATIC_VARIATIONS}", self.set_shifts())
                .replace("{COMMITHASH}", '"{}"'.format(self.commit_hash))
                .replace("{SETUP_IS_CLEAN}", self.setup_is_clean)
            )
        log.info("Code written to {}".format(self.executable))
        log.info("------------------------------------")
        log.info("Code Generation Report")
        log.info("------------------------------------")
        log.info("  Output path: {}".format(self.executable))
        log.info("  Total Number of Defines: {} ".format(self.number_of_defines))
        log.info("  Total Number of Outputs: {} ".format(self.number_of_outputs))
        log.info(
            "  Total Number of Output files: {} ".format(
                len(self._outputfiles_generated.keys())
            )
        )
        log.info("------------------------------------")

    def generate_main_code(self) -> Tuple[str, str]:
        """
        Generate the call commands for all the subsets. Additionally, generate all include statements for the main executable.
        Args:
            None
        Returns:
            Tuple, the generated calls and the generated includes
        """
        main_calls = ""
        for scope in self.scopes:
            main_calls += "    // {}\n".format(scope)
            main_calls += "".join(self.subset_calls[scope])
        main_includes = "".join(self.subset_includes)
        return main_calls, main_includes

    def get_cmake_path(self) -> str:
        """
        Get the path to the cmake file
        Args:
            None
        Returns:
            the path to the cmake file
        """
        return os.path.join(
            self.executable_name + "_generated_code", self.executable_name + ".cxx"
        )

    def generate_subsets(self, scope: str) -> None:
        """
        Generate the subsets for the given scope
        Args:
            scope: the scope to generate the subsets for
        Returns:
            None
        """
        log.debug(
            "Generating subsets for {} in scope {}".format(self.executable_name, scope)
        )
        log.debug(
            "Output folder: {}".format(
                os.path.join(self.output_folder, self.executable_name)
            )
        )
        log.debug("Producers: {}".format(self.configuration.producers[scope]))
        # in order to map the dfs correctly, we have to count the number of subset calls
        is_first = True
        counter = 0
        for producer in self.configuration.producers[scope]:
            subset = CodeSubset(
                file_name=producer.name,
                template=self.subset_template,
                producer=producer,
                scope=scope,
                folder=os.path.join(
                    self.output_folder, self.executable_name + "_generated_code"
                ),
                configuration_parameters=self.configuration.config_parameters[scope],
            )
            subset.create()
            subset.write()
            self.number_of_defines += subset.count
            log.debug(
                "Adding {} defines for {} in scope {}".format(
                    subset.count, producer.name, scope
                )
            )
            # two special cases:
            # 1. global scope: there we have to use df0 as the input df
            # 2. first call of all other scopes: we have to use the last global df as the input df
            if scope == self.global_scope and is_first:
                self.subset_calls[scope].append(
                    subset.call(inputscope="df0", outputscope=f"df{counter+1}_{scope}")
                )
            elif is_first:
                self.subset_calls[scope].append(
                    subset.call(
                        inputscope=f"df{self.main_counter[self.global_scope]}_{self.global_scope}",
                        outputscope=f"df{counter+1}_{scope}",
                    )
                )
            else:
                self.subset_calls[scope].append(
                    subset.call(
                        inputscope=f"df{counter}_{scope}",
                        outputscope=f"df{counter+1}_{scope}",
                    )
                )
            self.subset_includes.append(subset.include())
            self.main_counter[scope] += 1
            counter += 1
            is_first = False

    def generate_run_commands(self) -> str:
        """
        generate the dataframe snapshot commands for the main executable. A seperate output file is generated for each scope, that contains at least one output quantity.
        The process tracking is also generated here.

        Args:
            None
        Returns:
            str - the generated run commands

        """
        log.debug("Generating run commands")
        runcommands = ""
        for scope in self.scopes:
            outputset: List[str] = []
            for output in sorted(self.outputs[scope]):
                self.output_commands[scope].extend(output.get_leaves_of_scope(scope))
            if len(self.output_commands[scope]) > 0 and scope != self.global_scope:
                # if no output is produced by the scope, we do not create a corresponding output file
                self._outputfiles_generated[scope] = "outputpath_{scope}".format(
                    scope=scope
                )
                # convert output lists to a set to remove duplicates
                outputset = list(
                    set(
                        self.output_commands[scope]
                        + self.output_commands[self.global_scope]
                    )
                )
                # sort the output list to get alphabetical order of the output names
                outputset.sort()
                outputstring = '", "'.join(outputset)

                self.number_of_outputs += len(self.output_commands[scope])
                runcommands += "    auto {scope}_cutReport = df{counter}_{scope}.Report();\n".format(
                    scope=scope, counter=self.main_counter[scope]
                )
                runcommands += '    std::string {outputname} = std::regex_replace(std::string(output_path), std::regex("\\\\.root"), "_{scope}.root");\n'.format(
                    scope=scope, outputname=self._outputfiles_generated[scope]
                )
                runcommands += '    auto {scope}_result = df{counter}_{scope}.Snapshot("ntuple", {outputname}, {{"{outputstring}"}}, dfconfig);\n'.format(
                    scope=scope,
                    counter=self.main_counter[scope],
                    outputname=self._outputfiles_generated[scope],
                    outputstring=outputstring,
                )
        # add code for tracking the progress
        runcommands += self.set_process_tracking()
        # add code for the time taken for the dataframe setup
        runcommands += self.set_setup_printout()
        # add trigger of dataframe execution, for nonempty scopes
        for scope in self.scopes:
            if len(self.output_commands[scope]) > 0 and scope != self.global_scope:
                runcommands += f"    {scope}_result.GetValue();\n"
                runcommands += f'    Logger::get("main")->info("{scope}:");\n'
                runcommands += f"    {scope}_cutReport->Print();\n"
        log.info(
            "Output files generated for scopes: {}".format(
                self._outputfiles_generated.keys()
            )
        )

        return runcommands

    def set_debug_flag(self) -> str:
        """
        Set the debug flag in the template if the debug variable is set to true

        Returns:
            None
        """
        if self.debug:
            return "bool debug = true;"
        else:
            return "bool debug = false;"

    def set_shifts(self) -> str:
        """
        Set the shifts in the template if the debug variable is set to true

        Returns:
            None
        """
        shifts = "{"
        for scope in self._outputfiles_generated.keys():
            shifts += '{{ {outputname}, {{"'.format(
                outputname=self._outputfiles_generated[scope]
            )
            shifts += '", "'.join([s for s in self.configuration.shifts[scope]])
            shifts += '"} },'
        shifts = shifts[:-1] + "}"
        return shifts

    def set_output_quantities(self) -> str:
        """
        Set the output quantities in the template if the debug variable is set to true

        Returns:
            None
        """
        output_quantities = "{"
        for scope in self._outputfiles_generated.keys():
            output_quantities += '{{ {outputname}, {{"'.format(
                outputname=self._outputfiles_generated[scope]
            )
            output_quantities += '", "'.join(
                [q.name for q in self.configuration.outputs[scope]]
            )
            output_quantities += '"} },'
        output_quantities = output_quantities[:-1] + "}"
        return output_quantities

    def set_thead_flag(self, threads: int) -> None:
        """
        Set the multithreading flag in the template if the number of threads is greater than 1.

        Args:
            threads: The number of threads to be used.

        Returns:
            None
        """
        self.threads = threads

    def set_setup_printout(self) -> str:
        """
        adds the code for the timing information on the dataframe setup to the run commands.
        """
        printout = ""
        printout += '   Logger::get("main")->info("Finished Setup");\n'
        printout += '   Logger::get("main")->info("Runtime for setup (real time: {0:.2f}, CPU time: {1:.2f})",\n'
        printout += "                           timer.RealTime(), timer.CpuTime());\n"
        printout += "   timer.Continue();\n"
        printout += '   Logger::get("main")->info("Starting Evaluation");\n'

        return printout

    def set_process_tracking(self) -> str:
        """This function replaces the template placeholder for the process tracking with the correct process tracking.

        Returns:
            The code to be added to the template
        """
        tracking = ""
        scope = self.scopes[-1]
        tracking += "    ULong64_t {scope}_processed = 0;\n".format(scope=scope)
        tracking += "    std::mutex {scope}_bar_mutex;\n".format(scope=scope)
        tracking += "    auto c_{scope} = df{counter}_{scope}.Count();\n".format(
            counter=self.main_counter[scope], scope=scope
        )
        tracking += "    c_{scope}.OnPartialResultSlot(quantile, [&{scope}_bar_mutex, &{scope}_processed, &quantile, &nevents](unsigned int /*slot*/, ULong64_t /*_c*/) {{".format(
            scope=scope
        )
        tracking += (
            "\n        std::lock_guard<std::mutex> lg({scope}_bar_mutex);\n".format(
                scope=scope
            )
        )
        tracking += "        {scope}_processed += quantile;\n".format(scope=scope)
        tracking += (
            "        float percentage = 100 * {scope}_processed / nevents;\n".format(
                scope=scope
            )
        )
        tracking += '        Logger::get("main")->info("{{0:c}} / {{1:c}} ({{2:.2f}} %) Events processed ...", {scope}_processed, nevents, percentage);\n'.format(
            scope=scope
        )
        tracking += "    });\n"
        return tracking
