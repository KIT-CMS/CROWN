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
                            {"df": "df%i" % df_count}
                            if df_scope_count == 0
                            else {"df": "%s_df%i" % (scope, df_scope_count)}
                        )
                    )
                    + ";\n"
                )
                df_scope_count += 1
    commandlist += "    auto df_final = mt_df%i;" % df_scope_count
    return t.replace("    // {CODE_GENERATION}", commandlist)
