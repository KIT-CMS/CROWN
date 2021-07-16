import logging
from git import Repo

log = logging.getLogger(__name__)


class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


def fill_template(t, config):
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
            commandlist += (
                "    auto df%i = " % (df_count + 1)
                + call.format_map(
                    SafeDict(
                        {"df": "df%i" % df_count, "vec_open": "{", "vec_close": "}"}
                    )
                )
                + ";\n"
            )
            df_count += 1
            log.debug("Adding call for {}".format(producer.name))
            log.debug("|---> {}".format(commandlist.split("\n")[-2]))
    commandlist += "    auto global_df_final = df%i;\n" % df_count
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
        runcommands += '    auto %s_result = %s_df_final.Snapshot("ntuple", std::string(output_path) + "test_%s.root", %s, dfconfig);\n' % (
            scope,
            scope,
            scope,
            '{"'
            + '", "'.join(
                [
                    '", "'.join(q.get_leaves_of_scope(scope))
                    for q in config["output"][scope]
                ]
            )
            + '"}',
        )
    for scope in config["output"]:
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
    plain_output_list = (
        '{"' + '", "'.join([q.name for q in config["output"][scope]]) + '"}'
    )
    shiftset = set()
    for q in config["output"][scope]:
        for shift in q.get_shifts(scope):
            shiftset.add(shift)
    shiftlist = '{"' + '", "'.join(shiftset) + '"}'
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
        .replace(
            "{METADATAFILENAME}", 'std::string(output_path) + "test_%s.root"' % scope
        )
        .replace("{OUTPUT_QUANTITIES}", plain_output_list)
        .replace("{SYSTEMATIC_VARIATIONS}", shiftlist)
        .replace("{COMMITHASH}", '"%s"' % current_commit)
        .replace("{CLEANSETUP}", setup_is_clean)
    )
