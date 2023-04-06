from os import path, makedirs
import importlib
from code_generation.code_generation import CodeGenerator
from code_generation.friend_trees import FriendTreeConfiguration
import inspect


def run(args):
    # the unittest is based on the tau analysis config
    analysis_name = "unittest"
    available_samples = [
        "ggh_htautau",
        "ggh_hbb",
        "vbf_htautau",
        "vbf_hbb",
        "rem_htautau",
        "rem_hbb",
        "embedding",
        "embedding_mc",
        "ttbar",
        "diboson",
        "dyjets",
        "wjets",
        "data",
        "electroweak_boson",
    ]
    available_eras = ["2018"]
    available_scopes = ["et", "mt", "tt", "em", "mm", "ee"]

    ## setup variables
    shifts = set([shift for shift in args.shifts])
    sample_group = args.sample
    era = args.era
    scopes = list(set([scope.lower() for scope in args.scopes]))

    ## load config
    configname = args.config
    config = importlib.import_module(
        f"analysis_configurations.{analysis_name}.{configname}"
    )
    ## Setting up executables
    for scope in scopes:
        args.logger.info(f"Scope: {scope}")
        args.logger.info(f"Generating Friend code for {sample_group}...")
        args.logger.info(f"Configuration used: {config}")
        args.logger.info(f"Era: {era}")
        args.logger.info(f"Shifts: {shifts}")
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
        args.logger.info(f"Executable: {executable}")

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
