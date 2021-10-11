import argparse
from shutil import copyfile
from os import path, makedirs
import importlib
import logging
import logging.handlers
import sys
from code_generation.code_generation import fill_template

sys.dont_write_bytecode = True

parser = argparse.ArgumentParser(description="Generate the C++ code for a given config")
parser.add_argument("--template", type=str, help="Path to the template")
parser.add_argument("--output", type=str, help="Path to the output directory")
parser.add_argument("--analysis", type=str, help="Name of the analysis config")
parser.add_argument(
    "--channels",
    type=str,
    nargs="+",
    help='Channels to be activated. To select all, choose "auto"',
)
parser.add_argument(
    "--shifts",
    type=str,
    nargs="*",
    help='Shifts to be activated. To select all, choose "auto"',
)
parser.add_argument(
    "--samples",
    type=str,
    nargs="+",
    help='Samples to be processed. To select all, choose "auto"',
)
parser.add_argument(
    "--eras",
    type=str,
    nargs="+",
    help='Eras to be processed. To select all, choose "auto"',
)
parser.add_argument("--debug", type=str, help='set debug mode for building"')
args = parser.parse_args()
# Executables for each era and per following processes:
available_samples = ["ggh", "vbf", "rem_htt", "emb", "tt", "vv", "dy", "wj", "data"]
available_eras = ["2016", "2017", "2018"]


if "auto" in args.samples:
    args.samples = available_samples
if "auto" in args.eras:
    args.eras = available_eras

executables = []
for era in args.eras:
    for sample_group in args.samples:
        analysisname = args.analysis
        ## setup logging
        if not path.exists("generation_logs"):
            makedirs("generation_logs")
        handler = logging.handlers.WatchedFileHandler(
            "generation_logs/generation_{sample}.log".format(sample=sample_group), "w"
        )
        terminal_handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
        terminal_handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
        root = logging.getLogger()
        root.setLevel("INFO")
        if args.debug != "false":
            root.setLevel("DEBUG")
        root.addHandler(handler)
        root.addHandler(terminal_handler)

        ### Setting up executable
        analysis = importlib.import_module("config." + analysisname)
        executable = f"{analysisname}_{sample_group}_{era}.cxx"
        root.info("Generating code for {}...".format(sample_group))
        root.info("Configuration used: {}".format(analysis))
        config = analysis.build_config(
            era,
            sample_group,
            args.channels,
            args.shifts,
            available_samples,
            available_eras,
        )
        # fill code template and write executable
        with open(args.template, "r") as template_file:
            template = template_file.read()
        template = fill_template(template, config)
        template = (
            template.replace("{ANALYSISTAG}", '"Analysis=%s"' % analysisname)
            .replace("{ERATAG}", '"Era=%s"' % era)
            .replace("{SAMPLETAG}", '"Samplegroup=%s"' % sample_group)
        )
        with open(path.join(args.output, executable), "w") as executable_file:
            executable_file.write(template)
        executables.append(executable)

with open(path.join(args.output, "files.txt"), "w") as f:
    for filename in executables:
        # copyfile(args.template, path.join(args.output, filename))
        f.write(filename + "\n")
