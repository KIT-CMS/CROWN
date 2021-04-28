# Dummy Python generation script, just to showcase the cmake integration
# This script just copies the template around and does not modify it

import argparse
from shutil import copyfile
from os import path
import importlib

from code_generation.code_generation import fill_template

parser = argparse.ArgumentParser(description="Generate the C++ code for a given config")
parser.add_argument("--template", type=str, help="Path to the template")
parser.add_argument("--output", type=str, help="Path to the output directory")
parser.add_argument("--analysis", type=str, help="Name of the analysis config")
parser.add_argument(
    "--channels",
    type=str,
    nargs="*",
    help='Channels to be activated. To select all, choose "auto"',
)
parser.add_argument(
    "--shifts", type=str, help='Shifts to be activated. To select all, choose "auto"'
)
parser.add_argument(
    "--samples", type=str, help='Samples to be processed. To select all, choose "auto"'
)

args = parser.parse_args()
# Executables for each era and per following processes:
# ggH
# vbf
# remaining htt and hww?
# embedded
# tt
# vv and single top
# DY
# WJets
# data
executables = []
sample_groups = ["emb"]
for sample_group in sample_groups:
    executable = f"analysis_{sample_group}.cxx"
    analysis = importlib.import_module("config." + args.analysis)
    config = analysis.build_config()
    # modify config according to args
    # ...
    # process shifts
    # ...
    # fill code template and write executable
    with open(args.template, "r") as template_file:
        template = template_file.read()
    template = fill_template(template, config)
    with open(executable, "w") as executable_file:
        executable_file.write(template)
    executables.append(executable)

with open(path.join(args.output, "files.txt"), "w") as f:
    for filename in executables:
        # copyfile(args.template, path.join(args.output, filename))
        f.write(filename + "\n")
