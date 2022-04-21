import argparse
from os import path, makedirs
import importlib
import logging
import logging.handlers
import sys
from code_generation.code_generation import CodeGenerator

sys.dont_write_bytecode = True


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(","))


parser = argparse.ArgumentParser(description="Generate the C++ code for a given config")
parser.add_argument("--template", type=str, help="Path to the template")
parser.add_argument("--subset-template", type=str, help="Path to the subset template")
parser.add_argument("--output", type=str, help="Path to the output directory")
parser.add_argument("--analysis", type=str, help="Name of the analysis config")
parser.add_argument(
    "--scopes",
    type=str,
    action=SplitArgs,
    help="Channels to be activated. To select multiple scopes, provide a comma separated list.",
)
parser.add_argument(
    "--shifts",
    type=str,
    action=SplitArgs,
    help='Shifts to be activated. To select multiple shifts, provide a comma separated list. To select all, choose "all"',
)
parser.add_argument(
    "--sample",
    type=str,
    help="Sample to be processed",
)
parser.add_argument(
    "--era",
    type=str,
    help="Era to be processed",
)
parser.add_argument("--threads", type=int, help="number of threads to be used")
parser.add_argument("--debug", type=str, help="set debug mode for building")
args = parser.parse_args()

available_samples = [
    "ggh",
    "vbf",
    "rem_htt",
    "emb",
    "emb_mc",
    "tt",
    "vv",
    "dy",
    "wj",
    "data",
]
available_eras = ["2016", "2017", "2018"]
available_scopes = ["et", "mt", "tt", "em", "ee", "mm"]

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
    "generation_logs/generation_{era}_{sample}.log".format(
        era=era, sample=sample_group
    ),
    "w",
)
handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
root.addHandler(handler)
## load analysis config
analysisname = args.analysis
analysis = importlib.import_module("config." + analysisname)
## Setting up executable
executable = f"{analysisname}_{sample_group}_{era}.cxx"
root.info("Generating code for {}...".format(sample_group))
root.info("Configuration used: {}".format(analysis))
root.info("Era: {}".format(era))
root.info("Shifts: {}".format(shifts))
config = analysis.build_config(
    era,
    sample_group,
    scopes,
    shifts,
    available_samples,
    available_eras,
    available_scopes,
)
# ## fill code template and write executable
# with open(args.template, "r") as template_file:
#     template = template_file.read()
# # generate the code using the analysis config and the template
# template = fill_template(template, config)
# # set analysis, era and sampletags
# template = set_tags(template, analysisname, era, sample_group)
# # if the number of threads is greater than one, add the threading flag in the code
# template = set_thead_flag(template, args.threads)
# # generate the code for the process tracking in the df
# template = set_process_tracking(template, scopes)
# # set debug flag if running in debug mode
# template = set_debug_flag(template, args.debug)
# with open(path.join(args.output, executable), "w") as executable_file:
#     executable_file.write(template)

# create a CodeGenerator object
generator = CodeGenerator(
    main_template_path=args.template,
    sub_template_path=args.subset_template,
    configuration=config,
    executable_name=f"{analysisname}_{sample_group}_{era}",
    analysisname=analysisname,
    output_folder=args.output,
)
if args.debug == "true":
    generator.debug = True
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
