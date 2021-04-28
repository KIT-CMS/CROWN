import code_generation.producer as p


def fill_template(t, config):
    commandlist = "" #string to be placed into code template
    df_count = 0 #enumerate dataframes

    #get commands of producers and append to the command list
    for producer in config["producers"]:
        for call in getattr(p, producer).writecalls(config):
            # print(call)
            commandlist += (
                "  auto df%i = " % (df_count + 1)
                + call.format_map({"df": "df%i" % df_count})
                + ";\n"
            )
            df_count += 1
    commandlist += "  auto df_final = df%i;" % df_count
    return t.replace("  // {CODE_GENERATION}", commandlist)
