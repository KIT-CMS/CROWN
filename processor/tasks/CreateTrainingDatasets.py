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
from law.target.collection import \
    flatten_collections, NestedSiblingFileCollection
from ast import literal_eval
import re

# Task to create config files for the necessary datashards
# One config is created for each process (like "ff" and "NMSSM_240_125_60")
# Shards are NOT shared between eras and decay channels
class CreateTrainingDataShardConfig(Task):
    # Define luigi parameters
    eras = luigi.ListParameter(description="List of run eras")
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
            "/".join([
                self.dir_template.format(
                    era=era, 
                    channel=self.channel
                ),
                self.file_template.format(
                    process=process
                )
            ])
            for era in self.eras
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
        #Create list of missing processes, process classes and output targets
        missing_processes = list(compress(processes, missing_confs_bools))
        missing_process_classes = list(compress(process_classes, missing_confs_bools))
        missing_outputs = [
            output for output in flatten_collections(self.output()) if not output.exists()
        ]
        # Create new temporary data directory
        prefix = self.temporary_local_path("")
        os.makedirs(prefix, exist_ok=True)
        run_loc = "sm-htt-analysis"
        files = [
            "/".join([
                prefix,
                self.dir_template.format(
                    era = era,
                    channel = self.channel
                ),
                self.file_template.format(
                    process=process
                )
            ])
            for era in self.eras
            for process in missing_processes
        ]

        # Loop over eras
        for era in self.eras:
            # Create output dir for each era
            os.makedirs("/".join([
                prefix, 
                "{}_{}".format(era, self.channel)
            ]), exist_ok=True)
            if (self.channel != "em"):
                friends = "${SVFit_Friends} ${HHKinFit_Friends} ${FF_Friends}"
            else:
                friends = "${SVFit_Friends} ${HHKinFit_Friends}"
            # Write the config files for each missing datashard
            self.run_command(
                command=[
                    "python",
                    "ml_datasets/write_datashard_config.py",
                    "--era {}".format(era),
                    "--channel {}".format(self.channel),
                    "--processes {}".format(" ".join(missing_processes)),
                    "--process-classes {}".format(" ".join(missing_process_classes)),
                    "--base-path $ARTUS_OUTPUTS",
                    "--friend-paths {}".format(friends),
                    "--database $KAPPA_DATABASE",
                    "--tree-path {}_nominal/ntuple".format(self.channel),
                    "--event-branch event",
                    "--training-weight-branch training_weight",
                    "--output-config-base {}".format(prefix)
                ],
                run_location=run_loc,
                sourcescripts=[
                    "{}/utils/setup_python.sh".format(run_loc),
                    "{}/utils/setup_samples.sh {}".format(run_loc, era)
                ]
            )
        # self.run_command(command=["ls", "-R"], run_location=prefix)
        # Copy locally created files to remote storage
        for file_remote, file_local in zip(missing_outputs, files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        # Remove temporary data directory
        rmtree(prefix)


# Task to create root shards for the NN training
# One shard is created for each process (like "ff" and "NMSSM_240_125_60")
# Shards are NOT shared between eras and decay channels  
class CreateTrainingDataShard(HTCondorWorkflow, law.LocalWorkflow):
    # Define luigi parameters
    eras = luigi.ListParameter(description="List of run eras")
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
                "era": era,
                "channel": self.channel,
                "process": process,
                "fold": fold,
            } 
            for era in self.eras
            for process in processes
            for fold in ["0","1"]
        ]
        assert branches, \
            "There are no valid branches for this set of parameters: \
            \n{}".format(self)
        return branches
        
    # Set prerequisites of this task:
    # All configs for the dataset shards have to be completed
    def requires(self):
        return self.workflow_requires()

    def workflow_requires(self):
        processes, process_classes = zip(*self.processes_and_classes)
        requirements = super(CreateTrainingDataShard, self).workflow_requires()
        requirements_args = {
            "eras": self.eras, 
            "channel": self.channel,
            "processes_and_classes": self.processes_and_classes
        }
        requirements["CreateTrainingDataShardConfig"] = \
            CreateTrainingDataShardConfig(**requirements_args)
        return requirements

    # Define output targets. Task is considerd complete if all targets are present.
    def output(self):
        files = [
            "/".join([
                self.dir_template.format(
                    era=self.branch_data["era"], 
                    channel=self.branch_data["channel"]
                ),
                self.files_template.format(
                    process=self.branch_data["process"],
                    fold=self.branch_data["fold"]
                )
            ])
        ]
        targets = self.remote_targets(files)
        for target in targets:
            target.parent.touch()
        return targets

    def run(self):
        prefix = self.temporary_local_path("")
        run_loc = "sm-htt-analysis"
        # Create data directory
        data_dir = "/".join([
            prefix,
            self.dir_template.format(
                era=self.branch_data["era"],
                channel=self.branch_data["channel"]
            )
        ])
        os.makedirs(data_dir, exist_ok=True)
        files = [
            "/".join([
                data_dir,
                self.files_template.format(
                    process=self.branch_data["process"],
                    fold=self.branch_data["fold"]
                )
            ])
        ]
        # Get list of all shard config targets
        allbranch_targets = self.input()["CreateTrainingDataShardConfig"].targets
        # Filter targets by name in each branch
        filefilter = "/".join([
            self.dir_template.format(
                era = self.branch_data["era"], 
                channel = self.branch_data["channel"]
            ),
            self.files_template_in.format(
                process = self.branch_data["process"]
            )
        ])
        datashard_config = [
            target for target in allbranch_targets if filefilter in target.path
        ][0]

        # Copy filtered config file into data directory
        local_config_path = "/".join([
            data_dir,
            "{process}_datashard_config.yaml".format(
                process=self.branch_data["process"]
            )
        ])
        datashard_config.copy_to_local(local_config_path)

        # Create datashard based on config file
        self.run_command(
            [
                "python",
                "ml_datasets/create_training_datashard.py",
                "--config {}".format(local_config_path),
                "--process {}".format(self.branch_data["process"]),
                "--fold {}".format(self.branch_data["fold"]),
                "--output-dir {}/".format(data_dir)
            ],
            run_location=run_loc
        )
        # Copy locally created files to remote storage
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        # Remove data directory
        rmtree(prefix)

# Task to create config files for the NN trainings
# One config file is created for each valid combination of:
# era, channel, mass and batch
# For this the shards are combined and reused as necessary
class CreateTrainingConfig(Task, law.LocalWorkflow):
    # Define luigi parameters
    eras = luigi.ListParameter(description="List of run eras")
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

    # Set other variables and templates used by this class.
    dir_template = "{era}_{channel}_{mass}_{batch}"
    dir_template_in = "{era}_{channel}"
    file_template = "training_config.yaml"
    files_template_shard = "{process}_datashard_fold{fold}.root"
    files_template_shard_config = "{process}_datashard_config.yaml"

    # Function to check for invalid mass-batch combinations
    def valid_batches(self, mass):
        mass = str(mass)
        if mass in ["240","280"]:
            max_batch = 3
        elif mass in ["320","360","400","450"]:
            max_batch = 4
        elif mass in ["500","550","600"]:
            max_batch = 5
        elif mass in ["700","800","heavier"]:
            max_batch = 6
        elif mass in ["900","1000"]:
            max_batch = 7
        else:
            raise Exception("Provided mass {} is not valid.".format(mass))
        valid_batches = [batch for batch in self.batch_nums if int(batch) <= max_batch]
        return valid_batches

    # Create map for the branches of this task
    def create_branch_map(self):
        branches = [
            {
                "era": era,
                "channel": self.channel
            } 
            for era in self.eras
        ]
        assert branches, \
            "There are no valid branches for this set of parameters: \
            \n{}".format(self)
        return branches

    # Set prerequisites of this task:
    # All configs for the dataset shards have to be completed
    # All dataset shards have to be completed
    def requires(self):
        return self.workflow_requires()

    def workflow_requires(self):
        requirements = super(CreateTrainingConfig, self).workflow_requires()

        requirements_args = {
            "eras": self.eras, 
            "channel": self.channel,
            "processes_and_classes": self.processes_and_classes,
        }
        requirements["CreateTrainingDataShard"] = \
            CreateTrainingDataShard(**requirements_args)
        requirements["CreateTrainingDataShardConfig"] = \
            CreateTrainingDataShardConfig(**requirements_args)

        return requirements

    # Define output targets. Task is considerd complete if all targets are present.
    def output(self):
        files = [
            "/".join([
                self.dir_template.format(
                    era=self.branch_data["era"], 
                    channel=self.branch_data["channel"],
                    mass=mass,
                    batch=batch_num
                ),
                self.file_template
            ])
            for mass in self.masses
            for batch_num in self.valid_batches(mass)
        ]
        targets = self.remote_targets(files)
        for target in targets:
            target.parent.touch()
        return targets

    def run(self):
        processes, process_classes = zip(*self.processes_and_classes)
        prefix = self.temporary_local_path("")
        run_loc = "sm-htt-analysis"
        # Create data directory for each valid branch and mass/batch combination
        for mass in self.masses:
            for batch_num in self.valid_batches(mass):
                os.makedirs("/".join([
                    prefix, 
                    self.dir_template.format(
                        era=self.branch_data["era"], 
                        channel=self.branch_data["channel"],
                        mass=mass,
                        batch=batch_num
                    )
                ]), exist_ok=True)
                    
        files = [
            "/".join([
                prefix,
                self.dir_template.format(
                    era=self.branch_data["era"], 
                    channel=self.branch_data["channel"],
                    mass=mass,
                    batch=batch_num
                ),
                self.file_template
            ])
            for mass in self.masses
            for batch_num in self.valid_batches(mass)
        ]

        # Get list of all shard targets
        allbranch_shards = flatten_collections(
            self.input()["CreateTrainingDataShard"]["collection"]
        )
        # Filter shard targets by name in each branch and process
        filefilter_shard_strings = [
            "/".join([
                ".*",
                self.dir_template_in.format(
                    era = self.branch_data["era"], 
                    channel = self.branch_data["channel"]
                ),
                self.files_template_shard.format(
                    process = process,
                    fold = "."
                )
            ])
            for process in processes
        ]
        filefilter_shards = [re.compile(string) for string in filefilter_shard_strings]
        shards = [
            target 
            for target in allbranch_shards 
            for filt in filefilter_shards if filt.match(target.path)
        ]
        
        # Get list of all shard config targets
        allbranch_shardconfigs = self.input()["CreateTrainingDataShardConfig"].targets
        # Filter shard config targets by name in each branch and process
        filefilter_shard_config_strings = [
            "/".join([
                ".*",
                self.dir_template_in.format(
                    era = self.branch_data["era"], 
                    channel = self.branch_data["channel"]
                ),
                self.files_template_shard_config.format(
                    process = process
                )
            ])
            for process in processes
        ]
        filefilter_shard_configs = [
            re.compile(string) for string in filefilter_shard_config_strings
        ]
        shardconfigs = [
            target 
            for target in allbranch_shardconfigs 
            for filt in filefilter_shard_configs if filt.match(target.path)
        ]
        
        # Copy filtered shard and shard config files into data directory
        for copy_file in shardconfigs + shards:
            copy_file_name = os.path.basename(copy_file.path)
            copy_file_name_path = os.path.basename(os.path.dirname(copy_file.path))
            copy_file.copy_to_local("/".join([
                prefix,
                copy_file_name_path,
                copy_file_name
            ]))
        # Write training config file 
        self.run_command(
            [
                "python",
                "ml_datasets/write_training_config.py",
                "--era {}".format(self.branch_data["era"]),
                "--channel {}".format(self.channel),
                "--masses {}".format(" ".join(self.masses)),
                "--batches {}".format(" ".join(self.batch_nums)),
                "--dataset-dir {}/{}_{}".format(prefix, self.branch_data["era"], self.channel),
                "--processes {}".format(" ".join(processes)),
                "--training-template datasets/templates/{}_{}_training.yaml".format(
                    self.branch_data["era"], self.channel
                ),
                "--write-weights True",
                "--overwrite-configs {}".format("")
            ],
            run_location=run_loc,
            sourcescripts="{}/utils/setup_python.sh".format(run_loc)
        )
        # Copy locally created files to remote storage
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        # Remove data directories
        rmtree(prefix)

# Task to create config files for the NN trainings if all eras are used at once
# One config file is created for each valid combination of:
# channel, mass and batch
# For this the training configs of the individual eras are merged
class CreateTrainingConfigAllEras(Task, law.LocalWorkflow):
    # Define luigi parameters
    eras = luigi.ListParameter(description="List of run eras")
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

    # Set other variables and templates used by this class
    all_eras = ["2016", "2017", "2018"]
    dir_template = "{era}_{channel}_{mass}_{batch}"
    file_template = "training_config.yaml"

    # Function to check for invalid mass-batch combinations
    def valid_batches(self, mass):
        mass = str(mass)
        if mass in ["240","280"]:
            max_batch = 3
        elif mass in ["320","360","400","450"]:
            max_batch = 4
        elif mass in ["500","550","600"]:
            max_batch = 5
        elif mass in ["700","800","heavier"]:
            max_batch = 6
        elif mass in ["900","1000"]:
            max_batch = 7
        else:
            raise Exception("Provided mass {} is not valid.".format(mass))
        valid_batches = [batch for batch in self.batch_nums if int(batch) <= max_batch]
        return valid_batches

    # Create map for the branches of this task
    def create_branch_map(self):
        branches = [
            {
                "era": era,
                "mass": mass,
                "batch_num": batch_num
            } 
            for era in self.eras
            for mass in self.masses
            for batch_num in self.valid_batches(mass)
        ]
        assert branches, \
            "There are no valid branches for this set of parameters: \
            \n{}".format(self)
        return branches

    # Set prerequisites of this task:
    # Training configs for each era have to be completed
    def requires(self):
        return self.workflow_requires()

    def workflow_requires(self):
        requirements = super(CreateTrainingConfigAllEras, self).workflow_requires()
        requirements_args = {
            "eras": self.all_eras, 
            "channel": self.channel,
            "processes_and_classes": self.processes_and_classes,
            "masses": self.masses, 
            "batch_nums": self.batch_nums,
        }
        requirements["CreateTrainingConfig"] = CreateTrainingConfig(**requirements_args)
        return requirements

    # Define output targets. Task is considerd complete if all targets are present.
    def output(self):
        files = [
            "/".join([
                self.dir_template.format(
                    era=self.branch_data["era"], 
                    channel=self.channel,
                    mass=self.branch_data["mass"],
                    batch=self.branch_data["batch_num"]
                ),
                self.file_template
            ])
        ]
        targets = self.remote_targets(files)
        for target in targets:
            target.parent.touch()
        return targets
        
    def run(self):
        prefix = self.temporary_local_path("")
        run_loc = "sm-htt-analysis"
        # Create data directory
        data_dir = "/".join([
            prefix,
            self.dir_template.format(
                era=self.branch_data["era"],
                channel=self.channel,
                mass=self.branch_data["mass"],
                batch=self.branch_data["batch_num"]
            )
        ])
        os.makedirs("/".join([data_dir, "all_eras"]), exist_ok=True)
        files = [
            "/".join([
                data_dir,
                "all_eras",
                self.file_template
            ])
        ]
        # Get list of all training config targets
        allbranch_trainingconfigs = flatten_collections(
            self.input()["CreateTrainingConfig"]["collection"]
        )
        # Filter training config targets by name for every era
        filefilter_training_config_string = "/".join([
            ".*",
            self.dir_template.format(
                era = ".*", 
                channel = self.channel,
                mass = self.branch_data["mass"],
                batch = self.branch_data["batch_num"]
            ),
            self.file_template
        ])
        filefilter_training_config = re.compile(filefilter_training_config_string)
        trainingconfigs = [
            target 
            for target in allbranch_trainingconfigs 
            if filefilter_training_config.match(target.path)
        ]
        training_config_names = [os.path.basename(file.path) for file in trainingconfigs]

        # Copy filtered training config files into data directories
        for config, era in zip(trainingconfigs, self.all_eras):
            config_name = os.path.basename(config.path)
            config.copy_to_local(
                "/".join([data_dir, era, config_name])
            )
        # Combine the configs of each era
        self.run_command(
            [
                "python",
                "ml_datasets/create_combined_config.py",
                "--input-base-path {}".format(data_dir),
                "--output-dir {}".format(data_dir)
            ],
            run_location=run_loc
        )
        # Copy locally created files to remote storage
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        # Remove data directories
        rmtree(prefix)