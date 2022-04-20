from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
from typing import Any, Dict, List, Set, Union, Tuple
import os
import filecmp
from git import Repo, InvalidGitRepositoryError

from code_generation.producer import SafeDict, Producer, ProducerGroup

from code_generation.configuration import Configuration

log = logging.getLogger(__name__)


class CodeSubset(object):
    def __init__(
        self,
        filename: str,
        template: str,
        producer: Union[Producer, ProducerGroup],
        scope: str,
        folder: str,
        parameters: Dict[str, Any],
    ):
        self.filename = filename
        self.template = template
        self.producer = producer
        self.scope = scope
        self.name = producer.name + "_" + scope
        self.parameters = parameters
        self.count = 0
        self.folder = folder
        self.commands: List[str] = []
        self.headerfile = os.path.join(
            self.folder, "include", self.scope, "{}.hxx".format(self.filename)
        )
        self.sourcefile = os.path.join(
            self.folder, "src", self.scope, "{}.cxx".format(self.filename)
        )

    def create(self):
        """
        Create the code subset

        Args:
            None

        Returns:
            None
        """
        log.info("Creating code subset {}".format(self.name))
        log.info("Producer: {}".format(self.producer.name))
        log.info("Scope: {}".format(self.scope))
        self.producer.reserve_output(self.scope)
        # create the function calls for the producer
        for call in self.producer.writecalls(self.parameters, self.scope):
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
        Write the code subset to a file, both the header and the source. Before writing the files, check if they already exists, and if they exist and are not different, skip writing them. This is to avoid unnecessary recompilation, since the compiler will check the timestamps of the files.

        Args:
            None

        Returns:
            None
        """
        log.debug("Writing code subset {}".format(self.name))
        log.info("folder: {}, filename: {}".format(self.folder, self.filename))
        # write the header file if it does not exist or is different
        with open(self.headerfile + ".new", "w") as f:
            f.write(f"ROOT::RDF::RNode {self.name}(ROOT::RDF::RNode df);")
        if os.path.isfile(self.headerfile):
            if filecmp.cmp(self.headerfile + ".new", self.headerfile):
                log.info("--> Identical header file, skipping")
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
                log.info("--> Identical source file, skipping")
            else:
                os.rename(self.sourcefile + ".new", self.sourcefile)
        else:
            os.rename(self.sourcefile + ".new", self.sourcefile)

    def call(self, main_counter: int, scope: str, is_first: bool) -> str:
        """
        Return the call to the code subset

        Args:
            None

        Returns:
            str: the call to the code subset
        """
        if is_first:
            call = f"    auto df{main_counter+1}_{scope} = {self.name}(df0); \n"
        else:
            call = f"    auto df{main_counter+1}_{scope} = {self.name}(df{main_counter}_{scope}); \n"
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
    Class used to generate code from a given Configuration
    """

    def __init__(
        self,
        main_template_path: str,
        sub_template_path: str,
        configuration: Configuration,
        name: str,
        output_folder: str,
    ):
        self.main_template = self.load_template(main_template_path)
        self.subset_template = self.load_template(sub_template_path)
        self.configuration = configuration
        self.scopes = self.configuration.scopes
        self.outputs = self.configuration.outputs
        self.global_scope = self.configuration.global_scope
        self.name = name
        self.output_folder = output_folder
        self.executable = os.path.join(
            output_folder, self.name + "_generated_code", self.name + ".cxx"
        )
        self.metadata = ""
        self.subset_includes = []
        self.output_commands = {}
        self.subset_calls = {}
        self.main_counter = {}
        for scope in self.scopes:
            self.main_counter[scope] = 0
            self.subset_calls[scope] = []
            self.output_commands[scope] = []

    def generate_code(self) -> None:
        """
        Generate the code from the configuration and create the subsets. Run through the whole configuration and create a subset for each producer within the configuration.

        Start with the global scope and then all other scopes. All generated code is stored in the folder self.name.

        return None
        """
        # start with the global scope
        for subfolder in ["src", "include"]:
            for scope in self.scopes:
                if not os.path.exists(
                    os.path.join(self.name + "_generated_code", subfolder, scope)
                ):
                    os.makedirs(
                        os.path.join(self.name + "_generated_code", subfolder, scope)
                    )
        # self.generate_subsets(self.global_scope)
        for scope in self.scopes:
            self.generate_subsets(scope)

        calls, includes = self.generate_main_code()
        run_commands = self.generate_run_commands()
        nruns = self.generate_nruns()

        self.write_code(calls, includes, "", run_commands)

    def load_template(self, template_path: str) -> str:
        """
        Load the template from the given path
        :param template_path: the path to the template
        :return: the template
        """
        with open(template_path, "r") as template_file:
            template = template_file.read()
        return template

    def write_code(
        self, calls: str, includes: str, metadata: str, run_commands: str
    ) -> None:
        """
        Write the code to the output folder

        Args:
            main_calls: the main calls
            includes: the includes
            metadata: the metadata
            run_commands: the run commands

        Returns:
            None
        """
        with open(self.executable, "w") as f:
            f.write(
                self.main_template.replace("    // {CODE_GENERATION}", calls)
                .replace("// {INCLUDES}", includes)
                .replace("    // {RUN_COMMANDS}", run_commands)
            )

    def generate_main_code(self) -> Tuple(str):
        """
        Generate the main code from the configuration
        :return: the generated code
        """
        main_calls = ""
        for scope in self.scopes:
            main_calls += "    // {}\n".format(scope)
            main_calls += "".join(self.subset_calls[scope])
            if scope == self.global_scope:
                main_calls += "\n"
                for subscope in [
                    scope for scope in self.scopes if scope is not self.global_scope
                ]:
                    main_calls += "    auto  df0_{} = df{}_global;\n".format(
                        subscope, self.main_counter[self.global_scope]
                    )
        main_includes = "".join(self.subset_includes)
        return main_calls, main_includes

    def get_cmake_path(self):
        """
        Get the path to the cmake file
        :return: the path to the cmake file
        """
        return os.path.join(self.name + "_generated_code", self.name + ".cxx")

    def generate_subsets(self, scope) -> None:
        """
        Generate the subsets for the given scope
        :param scope: the scope to generate the subsets for
        :return: None
        """
        log.info("Generating subsets for {} in scope {}".format(self.name, scope))
        log.info(
            "Output folder: {}".format(os.path.join(self.output_folder, self.name))
        )
        log.info("Producers: {}".format(self.configuration.producers[scope]))
        is_first = True
        for producer in self.configuration.producers[scope]:
            subset = CodeSubset(
                filename=producer.name,
                template=self.subset_template,
                producer=producer,
                scope=scope,
                folder=os.path.join(self.output_folder, self.name + "_generated_code"),
                parameters=self.configuration.config_parameters[scope],
            )
            subset.create()
            subset.write()
            if scope == self.global_scope and is_first:
                self.subset_calls[scope].append(
                    subset.call(self.main_counter[scope], scope, True)
                )
            else:
                self.subset_calls[scope].append(
                    subset.call(self.main_counter[scope], scope, False)
                )
            self.subset_includes.append(subset.include())
            self.main_counter[scope] += 1
            is_first = False

    def generate_run_commands(self) -> str:
        """
        Write the snapshot command to the main code
        :return: None
        """
        log.info("Generating run commands")
        runcommands = ""
        for scope in self.scopes:
            for output in sorted(self.outputs[scope]):
                self.output_commands[scope].extend(output.get_leaves_of_scope(scope))
            if len(self.output_commands[scope]) > 0:
                outputname = "outputpath_{scope}".format(scope=scope)
                outputstring = '", "'.join(self.output_commands[scope])
                runcommands += "    auto {scope}_cutReport = df{counter}_{scope}.Report();\n".format(
                    scope=scope, counter=self.main_counter[scope]
                )
                runcommands += '    std::string {outputname} = std::regex_replace(std::string(output_path), std::regex("\\\\.root"), "_{scope}.root");\n'.format(
                    scope=scope, outputname=outputname
                )
                runcommands += '    auto {scope}_result = df{counter}_{scope}.Snapshot("ntuple", {outputname}, {{"{outputstring}"}}, dfconfig);\n'.format(
                    scope=scope,
                    counter=self.main_counter[scope],
                    outputname=outputname,
                    outputstring=outputstring,
                )
        # add trigger of dataframe execution, for nonempty scopes
        for scope in self.scopes:
            if len(self.output_commands[scope]) > 0:
                runcommands += f"    {scope}_result.GetValue();\n"
                runcommands += f'    Logger::get("main")->info("{scope}:");\n'
                runcommands += f"    {scope}_cutReport->Print();\n"
        return runcommands

    def generate_nruns(self) -> str:
        """
        Write the nruns command to the main code
        :return: nruns command
        """
        # calls = ""
        # for scope in self.scopes:
        #     calls += " + df{counter}_{scope}.GetNRuns()".format(scope=scope, counter=self.main_counter[scope])
        # nruns = (
        # "("
        # + "+".join(["%s_df_final.GetNRuns()" % scope for scope in config["producers"]])
        # + ")/%f" % len(config["producers"]))
        # return nruns
        pass

    def generate_metadata(self) -> None:
        """
        Write the metadata to the main code
        :return: None
        """
        pass


def set_tags(template: str, analysisname: str, era: str, sample_group: str) -> str:
    """
    Function used to set the tags in the template.

    Args:
        template: The template to be modified.
        analysisname: The name of the analysis.
        era: The era of the analysis.
        sample_group: The sample group of the analysis.

    Returns:
        The modified template.
    """
    return (
        template.replace("{ANALYSISTAG}", '"Analysis=%s"' % analysisname)
        .replace("{ERATAG}", '"Era=%s"' % era)
        .replace("{SAMPLETAG}", '"Samplegroup=%s"' % sample_group)
    )


def set_thead_flag(template: str, threads: int) -> str:
    """
    Set the multithreading flag in the template if the number of threads is greater than 1.

    Args:
        template: The template to be modified.
        threads: The number of threads to be used.

    Returns:
        The modified template.
    """
    if threads > 1:
        return template.replace(
            "// {MULTITHREADING}", "ROOT::EnableImplicitMT({});".format(threads)
        )
    else:
        return template.replace("// {MULTITHREADING}", "")


def set_process_tracking(template: str, scopes: List[str]) -> str:
    """This function replaces the template placeholder for the process tracking with the correct process tracking.

    Args:
        template: The template to be modified.
        scopes: The list of scopes to be used.

    Returns:
        The modified template.
    """
    tracking = ""
    for scope in scopes:
        tracking += "    ULong64_t {scope}_processed = 0;\n".format(scope=scope)
        tracking += "    std::mutex {scope}_bar_mutex;\n".format(scope=scope)
        tracking += "    auto c_{scope} = {scope}_df_final.Count();\n".format(
            scope=scope
        )
        tracking += "    c_{scope}.OnPartialResultSlot(quantile, [&{scope}_bar_mutex, &{scope}_processed, &quantile](unsigned int /*slot*/, ULong64_t /*_c*/) {{".format(
            scope=scope
        )
        tracking += (
            "\n        std::lock_guard<std::mutex> lg({scope}_bar_mutex);\n".format(
                scope=scope
            )
        )
        tracking += "        {scope}_processed += quantile;\n".format(scope=scope)
        tracking += '        Logger::get("main - {scope} scope")->info("{{}} Events processed ...", {scope}_processed);\n'.format(
            scope=scope
        )
        tracking += "    });\n"
    return template.replace("{PROGRESS_CALLBACK}", tracking)


def set_debug_flag(template: str, debug: str) -> str:
    """
    Set the debug flag in the template if the debug variable is set to true

    Args:
        template: The template to be modified.
        debug: The number of threads to be used.

    Returns:
        The modified template.
    """
    if debug == "true":
        return template.replace("    // {DEBUGLEVEL}", "    bool debug = true;")
    else:
        return template.replace("    // {DEBUGLEVEL}", "    bool debug = false;")
