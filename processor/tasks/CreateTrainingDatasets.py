# coding: utf-8

"""
Collection of tasks used to create training datasets and config files
for the NN trainings of the NMSSM analysis 
"""

import os
import luigi
import law
from shutil import rmtree
from framework import Task, HTCondorWorkflow
from itertools import compress
from law.target.collection import flatten_collections, NestedSiblingFileCollection
from ast import literal_eval
import re

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
        self.publish_message("File copy out start.")
        for file_remote, file_local in zip(missing_outputs, files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        self.publish_message("File copy out end.")
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
        }
        requirements["CreateTrainingDataShardConfig"] = CreateTrainingDataShardConfig(
            **requirements_args
        )
        return requirements

    def workflow_requires(self):
        processes, process_classes = zip(*self.processes_and_classes)
        requirements = {}
        requirements_args = {
            "era": self.era,
            "channel": self.channel,
            "processes_and_classes": self.processes_and_classes,
        }
        requirements["CreateTrainingDataShardConfig"] = CreateTrainingDataShardConfig(
            **requirements_args
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
        self.publish_message("File copy in start.")
        datashard_config.copy_to_local(local_config_path)
        self.publish_message("File copy in end.")

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
        self.publish_message("File copy out start.")
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        self.publish_message("File copy out end.")
        # Remove data directory
        rmtree(prefix)


# Task to create config files for the NN trainings
# One config file is created for each valid combination of:
# era, channel, mass and batch
# For this the shards are combined and reused as necessary
class CreateTrainingConfig(Task):
    # Define luigi parameters
    era = luigi.Parameter(description="Run era")
    channel = luigi.Parameter(description="Decay Channel")
    masses = luigi.ListParameter(
        description="List of mass hypotheses of heavy NMSSM Higgs boson."
    )
    batch_nums = luigi.ListParameter(
        description="List of groups of mass hypotheses of light NMSSM Higgs boson."
    )
    processes_and_classes = luigi.ListParameter(
        description="List of processes and their classes used in training"
    )

    all_eras = ["2016", "2017", "2018"]

    # Set other variables and templates used by this class.
    dir_template = "{era}_{channel}_{mass}_{batch}"
    dir_template_in = "{era}_{channel}"
    file_template = "training_config.yaml"
    files_template_shard = "{process}_datashard_fold{fold}.root"
    files_template_shard_config = "{process}_datashard_config.yaml"

    # Function to check for invalid mass-batch combinations
    def valid_batches(self, mass):
        mass = str(mass)
        if mass in ["240", "280"]:
            max_batch = 3
        elif mass in ["320", "360", "400", "450"]:
            max_batch = 4
        elif mass in ["500", "550", "600"]:
            max_batch = 5
        elif mass in ["700", "800", "heavier"]:
            max_batch = 6
        elif mass in ["900", "1000"]:
            max_batch = 7
        else:
            raise Exception("Provided mass {} is not valid.".format(mass))
        valid_batches = [batch for batch in self.batch_nums if int(batch) <= max_batch]
        return valid_batches

    # Set prerequisites of this task:
    # All configs for the dataset shards have to be completed
    # All dataset shards have to be completed
    def requires(self):
        if self.era == "all_eras":
            use_eras = self.all_eras
        else:
            use_eras = [self.era]

        requirements = {}
        for era in use_eras:
            requirements_args = {
                "era": era,
                "channel": self.channel,
                "processes_and_classes": self.processes_and_classes,
            }
            requirements[
                "CreateTrainingDataShard_{}".format(era)
            ] = CreateTrainingDataShard(**requirements_args)
            requirements[
                "CreateTrainingDataShardConfig_{}".format(era)
            ] = CreateTrainingDataShardConfig(**requirements_args)

        return requirements

    # Define output targets. Task is considerd complete if all targets are present.
    def output(self):
        files = [
            "/".join(
                [
                    self.dir_template.format(
                        era=self.era,
                        channel=self.channel,
                        mass=mass,
                        batch=batch_num,
                    ),
                    self.file_template,
                ]
            )
            for mass in self.masses
            for batch_num in self.valid_batches(mass)
        ]
        targets = self.remote_targets(files)
        for target in targets:
            target.parent.touch()
        return targets

    def run(self):
        if self.era == "all_eras":
            use_eras = self.all_eras
        else:
            use_eras = [self.era]
        processes, process_classes = zip(*self.processes_and_classes)
        prefix = self.temporary_local_path("")
        run_loc = "sm-htt-analysis"
        # Create data directory for each valid branch and mass/batch combination
        for era in use_eras + [self.era]:
            for mass in self.masses:
                for batch_num in self.valid_batches(mass):
                    os.makedirs(
                        "/".join(
                            [
                                prefix,
                                self.dir_template.format(
                                    era=era,
                                    channel=self.channel,
                                    mass=mass,
                                    batch=batch_num,
                                ),
                            ]
                        ),
                        exist_ok=True,
                    )

        files = [
            "/".join(
                [
                    prefix,
                    self.dir_template.format(
                        era=self.era,
                        channel=self.channel,
                        mass=mass,
                        batch=batch_num,
                    ),
                    self.file_template,
                ]
            )
            for mass in self.masses
            for batch_num in self.valid_batches(mass)
        ]

        # Get list of all shard targets
        allbranch_shards = {
            era: flatten_collections(
                self.input()["CreateTrainingDataShard_{}".format(era)]["collection"]
            )
            for era in use_eras
        }
        # Get prefix of remote storage for root shards
        remote_shard_base = {
            era: self.wlcg_path + os.path.dirname(allbranch_shards[era][0].path)
            for era in use_eras
        }

        # Get shard config targets
        branch_shardconfigs = {
            era: self.input()["CreateTrainingDataShardConfig_{}".format(era)].targets
            for era in use_eras
        }

        for era in use_eras:
            assert len(branch_shardconfigs[era]) == len(
                processes
            ), "There should be {} targets, but there are {}".format(
                len(processes),
                len(branch_shardconfigs[era]),
            )

        # Copy shard config files into data directory
        self.publish_message("File copy in start.")
        for era in use_eras:
            for copy_file in branch_shardconfigs[era]:
                copy_file_name = os.path.basename(copy_file.path)
                copy_file_name_path = os.path.basename(os.path.dirname(copy_file.path))
                copy_file.copy_to_local(
                    "/".join([prefix, copy_file_name_path, copy_file_name])
                )
        self.publish_message("File copy in end.")
        # Write training config file
        for era in use_eras:
            self.run_command(
                [
                    "python",
                    "ml_datasets/write_training_config.py",
                    "--era {}".format(era),
                    "--channel {}".format(self.channel),
                    "--masses {}".format(" ".join(self.masses)),
                    "--batches {}".format(" ".join(self.batch_nums)),
                    "--config-dir {}/{}_{}".format(prefix, era, self.channel),
                    "--dataset-dir {}".format(remote_shard_base[era]),
                    "--processes {}".format(" ".join(processes)),
                    "--training-template datasets/templates/{}_{}_training.yaml".format(
                        era, self.channel
                    ),
                    "--write-weights True",
                    "--overwrite-configs {}".format(""),
                ],
                run_location=run_loc,
                sourcescripts="{}/utils/setup_python.sh".format(run_loc),
            )
        if self.era == "all_eras":
            # Create data directory
            for mass in self.masses:
                for batch_num in self.valid_batches(mass):
                    conf_files = "/".join(
                        [
                            prefix,
                            self.dir_template.format(
                                era="{era}",
                                channel=self.channel,
                                mass=mass,
                                batch=batch_num,
                            ),
                        ]
                    )
                    print(conf_files)
                    # Combine the configs of each era
                    self.run_command(
                        [
                            "python",
                            "ml_datasets/create_combined_config.py",
                            "--input-path {}".format(conf_files),
                            "--output-dir {}".format(conf_files.format(era="all_eras")),
                        ],
                        run_location=run_loc,
                    )

        # Copy locally created files to remote storage
        self.publish_message("File copy out start.")
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        self.publish_message("File copy out end.")
        # Remove data directories
        rmtree(prefix)
