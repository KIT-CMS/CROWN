import argparse
import importlib
import os
import inspect
import logging
from code_generation.logger import CustomFormatter
import logging.handlers


class SplitArgs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values.split(","))


parser = argparse.ArgumentParser(description="Generate the C++ code for a given config")
parser.add_argument("--template", type=str, help="Path to the template")
parser.add_argument("--subset-template", type=str, help="Path to the subset template")
parser.add_argument("--output", type=str, help="Path to the output directory")
parser.add_argument("--analysis", type=str, help="Name of the analysis")
parser.add_argument("--config", type=str, help="Name of the config")
parser.add_argument(
    "--analysis-folder",
    type=str,
    help="The location of all analysis folders",
    default="analysis_configurations",
)
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
# the following arguments are only used for the production of friend trees
parser.add_argument("--friends", type=str, help="set production of friend trees")
parser.add_argument("--quantities-map", type=str, help="path to the quantities map")
args = parser.parse_args()

# find available analyses, every folder in analysis_configurations is an analysis
available_analysis = [
    os.path.basename(f.path)
    for f in os.scandir(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), args.analysis_folder)
    )
    if f.is_dir()
]

# sanitize arguments
available_analysis = [x.lower() for x in available_analysis]
analysis = args.analysis.lower()

# setup logging

## Setup Logging
root = logging.getLogger()
root.setLevel("INFO")
if args.debug == "true":
    root.setLevel("DEBUG")
## setup logging
if not os.path.exists("generation_logs"):
    os.makedirs("generation_logs")
terminal_handler = logging.StreamHandler()
fmt = "[%(levelname)-8s] [%(filename)-19s:%(lineno)-3d] %(message)s"
terminal_handler.setFormatter(CustomFormatter(fmt))
root.addHandler(terminal_handler)
handler = logging.handlers.WatchedFileHandler(
    f"generation_logs/generation_{args.era}_{args.sample}.log",
    "w",
)
handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
root.addHandler(handler)
args.logger = root
# first check if the analysis is available
if analysis not in available_analysis:
    raise ValueError(
        f"The analysis {analysis} is not available. Available analysiss are: {available_analysis}"
    )
configname = args.config
config = importlib.import_module(f"analysis_configurations.{analysis}.{configname}")
# check if the config is of type FriendTreeConfiguration
imported_members = [x[0] for x in inspect.getmembers(config, inspect.isclass)]
if (
    "FriendTreeConfiguration" in imported_members
    and "Configuration" in imported_members
):
    raise ValueError(
        f"Configuration {configname} contains both a Configuration and a FriendTreeConfiguration."
    )
elif "Configuration" in imported_members:
    basefile = "generate.py"
elif "FriendTreeConfiguration" in imported_members:
    if args.quantities_map is None:
        raise ValueError(
            "The quantity map is not provided. Please provide a path to the quantity map when using a FriendTreeConfiguration."
        )
    basefile = "generate_friends.py"
else:
    raise ValueError(
        f"Configuration {configname} does not contain a Configuration or a FriendTreeConfiguration."
    )

# load and run the generate.py function of the analysis
if not os.path.exists(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        args.analysis_folder,
        analysis,
        basefile,
    )
):
    raise ValueError(f"The generate.py file for analysis {analysis} does not exist")
generator = importlib.import_module(
    f"analysis_configurations.{analysis}.{basefile.replace('.py', '')}"
)
generator.run(args)
