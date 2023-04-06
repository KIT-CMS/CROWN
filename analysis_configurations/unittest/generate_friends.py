from os import path, makedirs
import importlib
import logging
import logging.handlers
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

    ## Setup Logging
    root = logging.getLogger()
    root.setLevel("INFO")
    if args.debug == "true":
        root.setLevel("DEBUG")
    ## setup logging
    if not path.exists("generation_logs"):
        makedirs("generation_logs")
    terminal_handler = logging.StreamHandler()
    terminal_handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    root.addHandler(terminal_handler)
    handler = logging.handlers.WatchedFileHandler(
        f"generation_logs/generation_{era}_{sample_group}.log",
        "w",
    )
    handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
    root.addHandler(handler)
    ## load config
    configname = args.config
    config = importlib.import_module(
        f"analysis_configurations.{analysis_name}.{configname}"
    )
    ## Setting up executables
    for scope in scopes:
        root.info(f"Scope: {scope}")
        root.info(f"Generating Friend code for {sample_group}...")
        root.info(f"Configuration used: {config}")
        root.info(f"Era: {era}")
        root.info(f"Shifts: {shifts}")
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
            analysis_name=f"{analysis_name}_{configname}",
            output_folder=args.output,
            threads=args.threads,
        )
        if args.debug == "true":
            generator.debug = True
        # generate the code
        generator.generate_code()

        executable = generator.get_cmake_path()
        logging.info(f"Executable: {executable}")

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
