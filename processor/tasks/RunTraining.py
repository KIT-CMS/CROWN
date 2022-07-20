# coding: utf-8

"""
Task to run NMSSM NN trainings
And Task to run all trainings for each mass,batch,channel combination
"""

import os
import luigi
import law
from shutil import rmtree
from rich.console import Console
from framework import Task, HTCondorWorkflow, PuppetMaster
from ast import literal_eval
from law.target.collection import flatten_collections
import re
from CreateTrainingDatasets import (
    CreateTrainingDataShardConfig,
    CreateTrainingDataShard,
)

law.contrib.load("tasks")  # to have the RunOnceTask class

try:
    current_width = os.get_terminal_size().columns
except OSError:
    current_width = 140
console = Console(width=current_width)


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
                "identifier": [self.channel],
            }
            requirements["CreateTrainingDataShard_{}".format(era)] = PuppetMaster(
                puppet_task=CreateTrainingDataShard(**requirements_args),
                identifier=[era, self.channel],
            )
            requirements["CreateTrainingDataShardConfig_{}".format(era)] = PuppetMaster(
                puppet_task=CreateTrainingDataShardConfig(**requirements_args),
                identifier=[era, self.channel],
            )
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
                self.requires()[
                    "CreateTrainingDataShard_{}".format(era)
                ].give_puppet_outputs()["collection"]
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
            era: self.requires()["CreateTrainingDataShardConfig_{}".format(era)]
            .give_puppet_outputs()
            .targets
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
        console.log("File copy in start.")
        for era in use_eras:
            for copy_file in branch_shardconfigs[era]:
                copy_file_name = os.path.basename(copy_file.path)
                copy_file_name_path = os.path.basename(os.path.dirname(copy_file.path))
                copy_file.copy_to_local(
                    "/".join([prefix, copy_file_name_path, copy_file_name])
                )
        console.log("File copy in end.")
        # Write training config file
        for era in use_eras:
            self.run_command(
                [
                    "python",
                    "ml_trainings/write_training_config.py",
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
                    # Combine the configs of each era
                    self.run_command(
                        [
                            "python",
                            "ml_trainings/create_combined_config.py",
                            "--input-path {}".format(conf_files),
                            "--output-dir {}".format(conf_files.format(era="all_eras")),
                        ],
                        run_location=run_loc,
                    )

        # Copy locally created files to remote storage
        console.log("File copy out start.")
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        console.log("File copy out end.")
        # Remove data directories
        rmtree(prefix)

# Task to run NN training (2 folds)
# One training is performed for each valid combination of:
# channel, mass, batch and fold
# The datashards are combined based on the training parameters
class RunTraining(HTCondorWorkflow, law.LocalWorkflow):
    # Define luigi parameters
    era = luigi.Parameter(description="Run era")
    channels = luigi.ListParameter(description="List of decay Channels")
    masses = luigi.ListParameter(
        description="List of mass hypotheses of heavy NMSSM Higgs boson."
    )
    batch_nums = luigi.ListParameter(
        description="List of groups of mass hypotheses of light NMSSM Higgs boson."
    )

    # Set other variables and templates used by this class.
    all_eras = ["2016", "2017", "2018"]
    dir_template = "{era}_{channel}_{mass}_{batch}"
    dir_template_in = "{era}_{channel}"
    file_template_shard = "{process}_datashard_fold{fold}.root"
    file_template_config = "training_config.yaml"
    file_templates = [
        "fold{fold}_keras_model.h5",
        "fold{fold}_keras_preprocessing.pickle",
        "fold{fold}_loss.pdf",
        "fold{fold}_loss.png",
    ]

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

    # Function to get relevant processes and their classes
    def get_process_tuple(self, channel, mass, batch_num):
        run_loc = "sm-htt-analysis"
        # Get processes and classes as list of tuples
        process_tuples = literal_eval(
            self.run_command(
                [
                    "python",
                    "utils/get_processes.py",
                    "--channel {}".format(channel),
                    "--mass {}".format(mass),
                    "--batch {}".format(batch_num),
                    "--training-z-estimation-method emb",
                    "--training-jetfakes-estimation-method ff",
                ],
                run_location=run_loc,
                sourcescripts="{}/utils/setup_python.sh".format(run_loc),
                collect_out=True,
                silent=True,
            )
        )
        return process_tuples

    # Create map for the branches of this task
    def create_branch_map(self):
        branches = [
            {"channel": channel, "mass": mass, "batch_num": batch_num, "fold": fold}
            for channel in self.channels
            for mass in self.masses
            for batch_num in self.valid_batches(mass)
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
    # All dataset shards have to be completed
    # All training config files have to be completed
    # The prerequisites are also dependant on whether all_eras is used
    def requires(self):
        channel = self.branch_data["channel"]
        mass = self.branch_data["mass"]
        batch_num = self.branch_data["batch_num"]
        requirements = {}

        processes_and_classes = self.get_process_tuple(
            channel,
            mass,
            batch_num,
        )
        requirements_conf = {
            "era": self.era,
            "channel": channel,
            "processes_and_classes": processes_and_classes,
            "masses": [mass],
            "batch_nums": [batch_num],
            "identifier": [channel],
        }
        requirements["CreateTrainingConfig_{}".format(channel)] = CreateTrainingConfig(
            **requirements_conf
        )
        if self.era == "all_eras":
            use_eras = self.all_eras
        else:
            use_eras = [self.era]
        # use_eras is set to be either given era, or a list of all eras
        for era in use_eras:
            requirements_shard = {
                "era": era,
                "channel": channel,
                "processes_and_classes": processes_and_classes,
                "identifier": [era, channel],
            }
            requirements[
                "CreateTrainingDataShard_{}_{}".format(
                    channel,
                    era,
                )
            ] = CreateTrainingDataShard(**requirements_shard)
        return requirements

    def workflow_requires(self):
        # Shortcut is necessary for Remote workflows without access to local filesystem
        if self.is_branch():
            return None
        requirements = super(RunTraining, self).workflow_requires()

        process_tuple = {
            channel: [
                self.get_process_tuple(channel, mass, batch_num)
                for mass in self.masses
                for batch_num in self.valid_batches(mass)
            ]
            for channel in self.channels
        }
        # Tasks of diffrent required channels are set up as seperate tasks, not branches
        # This is, because the process classes are dependant on the decay channel
        for channel in self.channels:
            # Get list of unique processes and their classes for each decay channel
            processes_and_classes = set(
                [
                    element
                    for element_list in process_tuple[channel]
                    for element in element_list
                ]
            )
            processes_and_classes = [list(elem) for elem in processes_and_classes]
            requirements_conf = {
                "era": self.era,
                "channel": channel,
                "processes_and_classes": processes_and_classes,
                "masses": self.masses,
                "batch_nums": self.batch_nums,
                "identifier": [channel],
            }
            # requires CreateTrainingConfig task if all_eras is not used
            requirements["CreateTrainingConfig_{}".format(channel)] = PuppetMaster(
                puppet_task=CreateTrainingConfig(**requirements_conf),
                identifier=[channel],
            )
            if self.era == "all_eras":
                use_eras = self.all_eras
            else:
                use_eras = [self.era]
            # use_eras is set to be either given era, or a list of all eras
            for era in use_eras:
                requirements_shard = {
                    "era": era,
                    "channel": channel,
                    "processes_and_classes": processes_and_classes,
                    "identifier": [era, channel],
                }
                requirements[
                    "CreateTrainingDataShard_{}_{}".format(channel, era)
                ] = PuppetMaster(
                    puppet_task=CreateTrainingDataShard(**requirements_shard),
                    identifier=[era, channel],
                )
        return requirements

    # Define output targets. Task is considerd complete if all targets are present.
    def output(self):
        files = [
            "/".join(
                [
                    self.dir_template.format(
                        era=self.era,
                        channel=self.branch_data["channel"],
                        mass=self.branch_data["mass"],
                        batch=self.branch_data["batch_num"],
                    ),
                    file_template.format(fold=self.branch_data["fold"]),
                ]
            )
            for file_template in self.file_templates
        ]
        targets = self.remote_targets(files)
        for target in targets:
            target.parent.touch()
        return targets

    def run(self):
        fold = self.branch_data["fold"]
        prefix = self.temporary_local_path("")
        run_loc = "sm-htt-analysis"
        # Create data directory for each branch
        local_dir = "/".join(
            [
                prefix,
                self.dir_template.format(
                    era=self.era,
                    channel=self.branch_data["channel"],
                    mass=self.branch_data["mass"],
                    batch=self.branch_data["batch_num"],
                ),
            ]
        )
        os.makedirs("/".join([local_dir, self.era]), exist_ok=True)
        if self.era == "all_eras":
            use_eras = self.all_eras
        else:
            use_eras = [self.era]
        files = [
            "/".join(
                [
                    local_dir,
                    self.era,
                    file_template.format(fold=self.branch_data["fold"]),
                ]
            )
            for file_template in self.file_templates
        ]
        # Get training config target
        branch_trainingconfigs = self.input()[
            "CreateTrainingConfig_{}".format(self.branch_data["channel"])
        ]
        assert (
            len(branch_trainingconfigs) == 1
        ), "There should be 1 target, but there are {}".format(
            len(branch_trainingconfigs)
        )
        trainingconfig = branch_trainingconfigs[0]

        # Get prefix of remote storage for root shards
        allbranch_shards = {
            era: flatten_collections(
                self.input()[
                    "CreateTrainingDataShard_{}_{}".format(
                        self.branch_data["channel"],
                        era,
                    )
                ]["collection"]
            )
            for era in use_eras
        }
        remote_shard_base = [
            self.wlcg_path
            + os.path.dirname(os.path.dirname(allbranch_shards[era][0].path))
            for era in use_eras
        ]
        assert (
            len(set(remote_shard_base)) == 1
        ), "Basepaths of different eras do not match"
        remote_shard_base = remote_shard_base[0]

        # Copy config to local
        trainingconfig_name = os.path.basename(trainingconfig.path)
        console.log("File copy in start.")
        trainingconfig.copy_to_local("/".join([local_dir, trainingconfig_name]))
        console.log("File copy in end.")

        # self.run_command(
        #     command=["ls", "-R"],
        #     run_location=local_dir
        # )
        # Set maximum number of threads (this number is somewhat inaccurate
        # as TensorFlow only abides by it for some aspects of the training)
        if "OMP_NUM_THREADS" in os.environ:
            max_threads = os.getenv("OMP_NUM_THREADS")
        else:
            max_threads = 12
        # Run NN training and save the model,
        # the preprocessing object and some images of the trasining process
        self.run_command(
            command=[
                "python",
                "ml_trainings/keras_training.py",
                "--data-dir {}".format(remote_shard_base),
                "--config-dir {}".format(local_dir),
                "--era-training {}".format(self.era),
                "--channel-training {}".format(self.branch_data["channel"]),
                "--fold {}".format(fold),
                "--balance-batches 1",
                "--max-threads {}".format(max_threads),
            ],
            run_location=run_loc,
        )
        # self.run_command(
        #     command=["ls", "-R"],
        #     run_location=prefix_w
        # )
        # Copy locally created files to remote storage
        console.log("File copy out start.")
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
        console.log("File copy out end.")
        # Remove data directories
        rmtree(prefix)


# Task to run all NN trainings for all decay channels,
# heavy higgs masses and groups of light higgs masses
class RunAllTrainings(Task, law.tasks.RunOnceTask):
    # Requires ALL trainings
    def requires(self):
        requirements = {
            "era": "all_eras",
            "channels": ["tt", "et", "mt"],
            "masses": [
                "240",
                "280",
                "320",
                "360",
                "400",
                "450",
                "500",
                "550",
                "600",
                "700",
                "800",
                "900",
                "1000",
                "heavier",
            ],
            "batch_nums": ["1", "2", "3", "4", "5", "6", "7"],
        }
        return PuppetMaster(puppet_task=RunTraining(**requirements))

    def run(self):
        console.log("All trainings are done!")
        self.mark_complete()
