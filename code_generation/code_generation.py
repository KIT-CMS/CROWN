import code_generation.producer as p


def fill_template(t, config):
    # shift producers
    for shift in config.keys():
        # shift names start with "_"
        if not shift.startswith("_"):
            continue
        for entry in config[shift]["shiftbase"]:
            getattr(p, entry).shift(shift)

    # generate list of commands
    commandlist = ""  # string to be placed into code template
    df_count = 0  # enumerate dataframes

    # get commands of producers and append to the command list
    for producer in config["producers"]:
        commandlist += "\n    //" + producer + "\n"
        for call in getattr(p, producer).writecalls(config):
            # print(call)
            commandlist += (
                "    auto df%i = " % (df_count + 1)
                + call.format_map(p.SafeDict({"df": "df%i" % df_count}))
                + ";\n"
            )
            df_count += 1
    commandlist += "    auto df_final = df%i;" % df_count
    return t.replace("    // {CODE_GENERATION}", commandlist)
