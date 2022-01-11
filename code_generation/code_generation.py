from __future__ import annotations  # needed for type annotations in > python 3.7

import logging
from typing import Any, Dict, List, Set

from git import Repo

from code_generation.producer import SafeDict

log = logging.getLogger(__name__)


def fill_template(t: str, config: Dict[Any, Any]) -> str:
    # generate list of commands
    commandlist = ""  # string to be placed into code template
    df_count = 0  # enumerate dataframes
    df_scope_count = 0  # enumerate dataframes in certain scopes
    # get commands of producers and append to the command list
    log.info("Generating commands ...")
    for producer in config["producers"]["global"]:
        producer.reserve_output("global")
        commandlist += "\n    //" + producer.name + "\n"
        for call in producer.writecalls(config, "global"):
            log.debug("Adding call for {}".format(producer.name))
            log.debug("Call: {}".format(call))
            commandlist += (
                "    auto df{} = ".format(df_count + 1)
                + call.format_map(
                    SafeDict(
                        {
                            "df": "df{}".format(df_count),
                            "vec_open": "{",
                            "vec_close": "}",
                        }
                    )
                )
                + ";\n"
            )
            df_count += 1
            log.debug("|---> {}".format(commandlist.split("\n")[-2]))
    commandlist += "    auto global_df_final = df{};\n".format(df_count)
    for scope in config["producers"]:
        if scope == "global":
            continue
        df_scope_count = 0
        for producer in config["producers"][scope]:
            producer.reserve_output(scope)
            commandlist += "\n    //" + producer.name + "\n"
            for call in producer.writecalls(config, scope):
                commandlist += (
                    "    auto %s_df%i = " % (scope, df_scope_count + 1)
                    + call.format_map(
                        SafeDict(
                            {"df": "global_df_final", "vec_open": "{", "vec_close": "}"}
                            if df_scope_count == 0
                            else {
                                "df": "%s_df%i" % (scope, df_scope_count),
                                "vec_open": "{",
                                "vec_close": "}",
                            }
                        )
                    )
                    + ";\n"
                )
                log.debug("Adding call for {}".format(producer.name))
                log.debug("|---> {}".format(commandlist.split("\n")[-2]))
                df_scope_count += 1
        commandlist += "    auto %s_df_final = %s_df%i;\n" % (
            scope,
            scope,
            df_scope_count,
        )
    commandlist += "\n"
    for scope in config["output"]:
        commandlist += "    auto %s_cutReport = %s_df_final.Report();\n" % (
            scope,
            scope,
        )
    runcommands = ""
    for scope in config["output"]:
        outputstring = (
            '{"'
            + '", "'.join(
                [
                    '", "'.join(q.get_leaves_of_scope(scope))
                    for q in config["output"][scope]
                ]
            )
            + '"}'
        )
        outputname = "outputpath_{scope}".format(scope=scope)
        runcommands += '    std::string {outputname} = std::regex_replace(std::string(output_path), std::regex("\\\\.root"), "_{scope}.root"); ;\n'.format(
            scope=scope, outputname=outputname
        )
        runcommands += '    auto {scope}_result = {scope}_df_final.Snapshot("ntuple", {outputname}, {outputstring}, dfconfig);\n'.format(
            scope=scope, outputname=outputname, outputstring=outputstring
        )
    for scope in config["output"]:
        runcommands += "{PROGRESS_CALLBACK}\n"
        runcommands += "    %s_result.GetValue();\n" % scope
        runcommands += '    Logger::get("main")->info("%s:");\n' % scope
        runcommands += "    %s_cutReport->Print();\n" % scope
    nruns = (
        "("
        + "+".join(["%s_df_final.GetNRuns()" % scope for scope in config["producers"]])
        + ")/%f" % len(config["producers"])
    )
    log.info("Finished generating code.")
    log.info("Prepare meta data.")
    plain_output_lists: List[str] = []
    for scope in config["output"]:
        outputname = "outputpath_{scope}".format(scope=scope)
        plain_output_lists.append(
            '{{ {outputname}, {{"'.format(outputname=outputname)
            + '", "'.join([q.name for q in config["output"][scope]])
            + '"} }'
        )
    output_lists = "{" + " , ".join(plain_output_lists) + " }"
    shiftset: Set[str] = set()
    shiftlists: List[str] = []
    for scope in config["output"]:
        outputname = "outputpath_{scope}".format(scope=scope)
        for q in config["output"][scope]:
            for shift in q.get_shifts(scope):
                shiftset.add(shift)
        shiftlists.append(
            '{{ {outputname}, {{"'.format(outputname=outputname)
            + '", "'.join(shiftset)
            + '"} }'
        )
    shifts = "{" + " , ".join(shiftlists) + " }"
    try:
        repo = Repo("../../CROWN")
        current_commit = repo.head.commit
        setup_is_clean = "false" if repo.is_dirty() else "true"
    except:
        current_commit = "undefined"
        setup_is_clean = "false"
    log.info("Finished preparing meta data.")
    return (
        t.replace("    // {CODE_GENERATION}", commandlist)
        .replace("    // {RUN_COMMANDS}", runcommands)
        .replace("{NRUNS}", nruns)
        .replace("{OUTPUT_QUANTITIES}", output_lists)
        .replace("{SYSTEMATIC_VARIATIONS}", shifts)
        .replace("{COMMITHASH}", '"%s"' % current_commit)
        .replace("{CLEANSETUP}", setup_is_clean)
    )


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


def set_process_tracking(template: str, channels: List[str]) -> str:
    """This function replaces the template placeholder for the process tracking with the correct process tracking.

    Args:
        template: The template to be modified.
        channels: The list of channels to be used.

    Returns:
        The modified template.
    """
    tracking = ""
    for channel in channels:
        tracking += "    ULong64_t {ch}_processed = 0;\n".format(ch=channel)
        tracking += "    std::mutex {ch}_bar_mutex;\n".format(ch=channel)
        tracking += "    auto c = {ch}_df_final.Count();\n".format(ch=channel)
        tracking += "    c.OnPartialResultSlot(quantile, [&{ch}_bar_mutex, &{ch}_processed, &quantile](unsigned int /*slot*/, ULong64_t /*_c*/) {{".format(
            ch=channel
        )
        tracking += (
            "\n        std::lock_guard<std::mutex> lg({ch}_bar_mutex);\n".format(
                ch=channel
            )
        )
        tracking += "        {ch}_processed += quantile;\n".format(ch=channel)
        tracking += '        Logger::get("main - {ch} Channel")->info("{{}} Events processed ...", {ch}_processed);\n'.format(
            ch=channel
        )
        tracking += "    });\n"
    return template.replace("{PROGRESS_CALLBACK}", tracking)
