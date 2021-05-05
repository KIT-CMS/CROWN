class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


def fill_template(t, config):
    # generate list of commands
    commandlist = ""  # string to be placed into code template
    df_count = 0  # enumerate dataframes
    df_scope_count = 0  # enumerate dataframes in certain scopes

    # get commands of producers and append to the command list
    for producer in config["producers"]["global"]:
        commandlist += "\n    //" + producer.name + "\n"
        for call in producer.writecalls(config, "global"):
            commandlist += (
                "    auto df%i = " % (df_count + 1)
                + call.format_map(SafeDict({"df": "df%i" % df_count}))
                + ";\n"
            )
            df_count += 1
    commandlist += "    auto global_df_final = df%i;\n" % df_count
    for scope in config["producers"]:
        if scope == "global":
            continue
        df_scope_count = 0
        for producer in config["producers"][scope]:
            commandlist += "\n    //" + producer.name + "\n"
            for call in producer.writecalls(config, scope):
                commandlist += (
                    "    auto %s_df%i = " % (scope, df_scope_count + 1)
                    + call.format_map(
                        SafeDict(
                            {"df": "global_df_final"}
                            if df_scope_count == 0
                            else {"df": "%s_df%i" % (scope, df_scope_count)}
                        )
                    )
                    + ";\n"
                )
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
                    '", "'.join(q.get_leafs_of_scope(scope))
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
    return (
        t.replace("    // {CODE_GENERATION}", commandlist)
        .replace("    // {RUN_COMMANDS}", runcommands)
        .replace("{NRUNS}", nruns)
    )
