class SafeDict(dict):
    def __missing__(self, key):
        return "{" + key + "}"


def fill_template(t, config):
    # generate list of commands
    commandlist = ""  # string to be placed into code template
    df_count = 0  # enumerate dataframes

    # get commands of producers and append to the command list
    for producer in config["producers"]:
        commandlist += "\n    //" + producer.name + "\n"
        for call in producer.writecalls(config):
            # print(call)
            commandlist += (
                "    auto df%i = " % (df_count + 1)
                + call.format_map(SafeDict({"df": "df%i" % df_count}))
                + ";\n"
            )
            df_count += 1
    commandlist += "    auto df_final = df%i;" % df_count
    return t.replace("    // {CODE_GENERATION}", commandlist)
