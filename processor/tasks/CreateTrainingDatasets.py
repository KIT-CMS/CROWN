# coding: utf-8

"""
Collection of tasks used to create training datasets and config files
for the NN trainings of the NMSSM analysis 
"""

import os
import luigi
import law
from shutil import rmtree
from rich.console import Console
from framework import Task, HTCondorWorkflow, PuppetMaster
from itertools import compress
from law.target.collection import flatten_collections, NestedSiblingFileCollection
from ast import literal_eval
import re

try:
    current_width = os.get_terminal_size().columns
except OSError:
    current_width = 140
console = Console(width=current_width)

# Task to create config files for the necessary datashards
# One config is created for each process (like "ff" and "NMSSM_240_125_60")
# Shards are NOT shared between eras and decay channels
class CreateTrainingDataShardConfig(Task):
    # Define luigi parameters
    era = luigi.Parameter(description="Run era")
    channel = luigi.Parameter(description="Decay Channel")
    processes_and_classes = luigi.ListParameter(
        description="List of processes and their classes used in training"
    )

    # Set other variables and templates used by this class.
    dir_template = "{era}_{channel}"
    file_template = "{process}_datashard_config.yaml"

    # Define output targets.
    # Task is considerd complete if all targets are present.
    def output(self):
        processes, process_classes = zip(*self.processes_and_classes)
        files = [
            "/".join(
                [
                    self.dir_template.format(era=self.era, channel=self.channel),
                    self.file_template.format(process=process),
                ]
            )
            for process in processes
        ]
        targets = self.remote_targets(files)
        for target in flatten_collections(targets):
            target.parent.touch()
        # Wrap targets into a nested sibling file collection.
        targets_collection = NestedSiblingFileCollection(targets)
        return targets_collection

    def run(self):
        processes, process_classes = zip(*self.processes_and_classes)
        # Check for each output target if it already exists
        missing_confs_bools = [
            not target.exists() for target in flatten_collections(self.output())
        ]
        # Create list of missing processes, process classes and output targets
        missing_processes = list(compress(processes, missing_confs_bools))
        missing_process_classes = list(compress(process_classes, missing_confs_bools))
        missing_outputs = list(
            compress(flatten_collections(self.output()), missing_confs_bools)
        )
        # Create new temporary data directory
        prefix = self.temporary_local_path("")
        os.makedirs(prefix, exist_ok=True)
        run_loc = "sm-htt-analysis"
        files = [
            "/".join(
                [
                    prefix,
                    self.dir_template.format(era=self.era, channel=self.channel),
                    self.file_template.format(process=process),
                ]
            )
            for process in missing_processes
        ]

        # Create output dir for each era
        os.makedirs(
            "/".join([prefix, "{}_{}".format(self.era, self.channel)]), exist_ok=True
        )
        if self.channel != "em":
            friends = "${SVFit_Friends} ${HHKinFit_Friends} ${FF_Friends}"
        else:
            friends = "${SVFit_Friends} ${HHKinFit_Friends}"
        # Write the config files for each missing datashard
        self.run_command(
            command=[
                "python",
                "ml_datasets/write_datashard_config.py",
                "--era {}".format(self.era),
                "--channel {}".format(self.channel),
                "--processes {}".format(" ".join(missing_processes)),
                "--process-classes {}".format(" ".join(missing_process_classes)),
                "--base-path $ARTUS_OUTPUTS",
                "--friend-paths {}".format(friends),
                "--database $KAPPA_DATABASE",
                "--tree-path {}_nominal/ntuple".format(self.channel),
                "--event-branch event",
                "--training-weight-branch training_weight",
                "--output-config-base {}".format(prefix),
            ],
            run_location=run_loc,
            sourcescripts=[
                "{}/utils/setup_python.sh".format(run_loc),
                "{}/utils/setup_samples.sh {}".format(run_loc, self.era),
            ],
        )
        # self.run_command(command=["ls", "-R"], run_location=prefix)
        # Copy locally created files to remote storage
        console.log("File copy out start.")
        for file_remote, file_local in zip(missing_outputs, files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        console.log("File copy out end.")
        # Remove temporary data directory
        rmtree(prefix)


# Task to create root shards for the NN training
# One shard is created for each process (like "ff" and "NMSSM_240_125_60")
# Shards are NOT shared between eras and decay channels
class CreateTrainingDataShard(HTCondorWorkflow, law.LocalWorkflow):
    # Define luigi parameters
    era = luigi.Parameter(description="Run era")
    channel = luigi.Parameter(description="Decay Channel")
    processes_and_classes = luigi.ListParameter(
        description="List of processes and their classes used in training"
    )

    # Set other variables and templates used by this class.
    dir_template = "{era}_{channel}"
    files_template = "{process}_datashard_fold{fold}.root"
    files_template_in = "{process}_datashard_config.yaml"

    # Create map for the branches of this task.
    # Each branch can be run as an individula job
    def create_branch_map(self):
        processes, process_classes = zip(*self.processes_and_classes)
        branches = [
            {
                "channel": self.channel,
                "process_and_class": process_and_class,
                "fold": fold,
            }
            for process_and_class in self.processes_and_classes
            for fold in ["0", "1"]
        ]
        assert (
            branches
        ), "There are no valid branches for this set of parameters: \
            \n{}".format(
            self
        )
        return branches

    # Set prerequisites of this task:
    # All configs for the dataset shards have to be completed
    def requires(self):
        processes, process_classes = zip(*self.processes_and_classes)
        requirements = {}
        requirements_args = {
            "era": self.era,
            "channel": self.channel,
            "processes_and_classes": [self.branch_data["process_and_class"]],
            "identifier": [self.era, self.channel],
        }
        requirements["CreateTrainingDataShardConfig"] = CreateTrainingDataShardConfig(
            **requirements_args
        )
        return requirements

    def workflow_requires(self):
        # Shortcut is necessary for Remote workflows without access to local filesystem
        if self.is_branch():
            return None
        processes, process_classes = zip(*self.processes_and_classes)
        requirements = {}
        requirements_args = {
            "era": self.era,
            "channel": self.channel,
            "processes_and_classes": self.processes_and_classes,
            "identifier": [self.era, self.channel],
        }
        requirements["CreateTrainingDataShardConfig"] = PuppetMaster(
            puppet_task=CreateTrainingDataShardConfig(**requirements_args),
            identifier=[self.era, self.channel],
        )
        return requirements

    # Define output targets. Task is considerd complete if all targets are present.
    def output(self):
        process, class_name = self.branch_data["process_and_class"]
        files = [
            "/".join(
                [
                    self.dir_template.format(
                        era=self.era, channel=self.branch_data["channel"]
                    ),
                    self.files_template.format(
                        process=process,
                        fold=self.branch_data["fold"],
                    ),
                ]
            )
        ]
        targets = self.remote_targets(files)
        for target in targets:
            target.parent.touch()
        return targets

    def run(self):
        process, class_name = self.branch_data["process_and_class"]
        prefix = self.temporary_local_path("")
        run_loc = "sm-htt-analysis"
        # Create data directory
        data_dir = "/".join(
            [
                prefix,
                self.dir_template.format(
                    era=self.era, channel=self.branch_data["channel"]
                ),
            ]
        )
        os.makedirs(data_dir, exist_ok=True)
        files = [
            "/".join(
                [
                    data_dir,
                    self.files_template.format(
                        process=process,
                        fold=self.branch_data["fold"],
                    ),
                ]
            )
        ]
        # Get config target
        allbranch_targets = self.input()["CreateTrainingDataShardConfig"].targets
        assert (
            len(allbranch_targets) == 1
        ), "There should be 1 target, but there are {}".format(len(allbranch_targets))
        datashard_config = allbranch_targets[0]

        # Copy config file into data directory
        local_config_path = "/".join(
            [
                data_dir,
                "{process}_datashard_config.yaml".format(process=process),
            ]
        )
        console.log("File copy in start.")
        datashard_config.copy_to_local(local_config_path)
        console.log("File copy in end.")

        # Create datashard based on config file
        # ROOT_XRD_QUERY_READV_PARAMS=0 is necessary for streaming root files from dCache
        self.run_command(
            [
                "ROOT_XRD_QUERY_READV_PARAMS=0",
                "python",
                "ml_datasets/create_training_datashard.py",
                "--config {}".format(local_config_path),
                "--process {}".format(process),
                "--fold {}".format(self.branch_data["fold"]),
                "--output-dir {}/".format(data_dir),
            ],
            run_location=run_loc,
        )
        # Copy locally created files to remote storage
        console.log("File copy out start.")
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        console.log("File copy out end.")
        # Remove data directory
        rmtree(prefix)
