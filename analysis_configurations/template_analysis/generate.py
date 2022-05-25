from os import path, makedirs
import importlib
import logging
import logging.handlers
from code_generation.code_generation import CodeGenerator


def run(args):

    analysis_name = "template_analysis"
    available_samples = [
        "data",
        "embedding",
        "ttbar",
        "dyjets",
        "wjets",
        "diboson",
    ]
    available_eras = ["2018"]
    available_scopes = ["mm"]

    ## setup variables
    shifts = set([shift.lower() for shift in args.shifts])
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
    ## Setting up executable
    executable = f"{configname}_{sample_group}_{era}.cxx"
    root.info(f"Generating code for {sample_group}...")
    root.info(f"Configuration used: {config}")
    root.info(f"Era: {era}")
    root.info(f"Shifts: {shifts}")
    config = config.build_config(
        era,
        sample_group,
        scopes,
        shifts,
        available_samples,
        available_eras,
        available_scopes,
    )
    # create a CodeGenerator object
    generator = CodeGenerator(
        main_template_path=args.template,
        sub_template_path=args.subset_template,
        configuration=config,
        executable_name=f"{configname}_{sample_group}_{era}",
        analysis_name=f"{analysis_name}_{configname}",
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
