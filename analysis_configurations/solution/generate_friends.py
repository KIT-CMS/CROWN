from os import path, makedirs
import importlib
from code_generation.code_generation import CodeGenerator
from code_generation.friend_trees import FriendTreeConfiguration
import inspect


def run(args):
    ### SET these variables according to your analysis ###
    analysis_name = "solution"
    available_samples = [
        "ggh_htautau",
        "ggh_hbb",
        "vbf_htautau",
        "vbf_hbb",
        "rem_htautau",
        "rem_hbb",
        "embedding",
        "embedding_mc",
        "singletop",
        "ttbar",
        "diboson",
        "dyjets",
        "wjets",
        "data",
        "electroweak_boson",
    ]
    available_eras = ["2016preVFP", "2016postVFP", "2017", "2018","2012"]
    available_scopes = ["et", "mt", "tt", "em", "ee", "mm"]
    ######################################################

    ## setup variables
    shifts = set([shift.lower() for shift in args.shifts])
    sample_group = args.sample
    era = args.era
    scopes = list(set([scope.lower() for scope in args.scopes]))

    ## load config
    configname = args.config
    config = importlib.import_module(
        f"analysis_configurations.{analysis_name}.{configname}"
    )
    # check if the config is of type FriendTreeConfiguration
    imported_members = [x[0] for x in inspect.getmembers(config, inspect.isclass)]
    if (
        "FriendTreeConfiguration" in imported_members
        and "Configuration" in imported_members
    ):
        raise ValueError(
            f"Configuration {configname} contains both a Configuration and a FriendTreeConfiguration."
        )
    elif "FriendTreeConfiguration" not in imported_members:
        raise ValueError(
            f"Configuration {configname} is not a FriendTreeConfiguration."
        )
    ## Setting up executable
    executable = f"{configname}_{sample_group}_{era}.cxx"
    args.logger.info(f"Generating code for {sample_group}...")
    args.logger.info(f"Configuration used: {config}")
    args.logger.info(f"Era: {era}")
    args.logger.info(f"Shifts: {shifts}")
    args.logger.info(f"Scopes: {scopes}")
    for scope in scopes:
        code_generation_config = config.build_config(
            era,
            sample_group,
            [scope],
            shifts,
            available_samples,
            available_eras,
            available_scopes,
            args.quantities_map,
        )
        # check if the config is of type FriendTreeConfiguration
        if not isinstance(code_generation_config, FriendTreeConfiguration):
            raise ValueError(
                f"Configuration {configname} is not a FriendTreeConfiguration."
            )
        # create a CodeGenerator object
        generator = CodeGenerator(
            main_template_path=args.template,
            sub_template_path=args.subset_template,
            configuration=code_generation_config,
            executable_name=f"{configname}_{sample_group}_{era}_{scope}",
            analysis_name=analysis_name,
            config_name=configname,
            output_folder=args.output,
            threads=args.threads,
        )
        if args.debug == "true":
            generator.debug = True
        # generate the code
        generator.generate_code()

        executable = generator.get_cmake_path()

        # append the executable name to the files.txt file
        # if the file does not exist, create it
        if not path.exists(path.join(args.output, "files.txt")):
            with open(path.join(args.output, "files.txt"), "w") as f:
                f.write(f"{executable}\n")
        else:
            with open(path.join(args.output, "files.txt"), "r+") as f:
                for line in f:
                    if executable in line:
                        break
                else:
                    f.write(f"{executable}\n")
