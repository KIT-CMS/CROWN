# coding: utf-8

"""
Task to run NMSSM NN trainings
And Task to run all trainings for each mass,batch,channel combination
"""

import os
import luigi
import law
from shutil import rmtree
from framework import Task, HTCondorWorkflow
from ast import literal_eval
from law.target.collection import flatten_collections
import re
from CreateTrainingDatasets import (
    CreateTrainingConfig,
    CreateTrainingConfigAllEras,
    CreateTrainingDataShard,
)

law.contrib.load("tasks")  # to have the RunOnceTask class

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
        return self.workflow_requires()

    def workflow_requires(self):
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
                "eras": [self.era],
                "channel": channel,
                "processes_and_classes": processes_and_classes,
                "masses": self.masses,
                "batch_nums": self.batch_nums,
            }
            if self.era == "all_eras":
                # requires CreateTrainingConfigAllEras task if all_eras is used
                requirements[
                    "CreateTrainingConfig_{}".format(channel)
                ] = CreateTrainingConfigAllEras(**requirements_conf)

                use_eras = self.all_eras
            else:
                # requires CreateTrainingConfig task if all_eras is not used
                requirements[
                    "CreateTrainingConfig_{}".format(channel)
                ] = CreateTrainingConfig(**requirements_conf)
                use_eras = [self.era]
            # use_eras is set to be either given era, or a list of all eras
            requirements_shard = {
                "eras": use_eras,
                "channel": channel,
                "processes_and_classes": processes_and_classes,
            }
            requirements[
                "CreateTrainingDataShard_{}".format(channel)
            ] = CreateTrainingDataShard(**requirements_shard)
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
        data_dir = "/".join(
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
        os.makedirs("/".join([data_dir, self.era]), exist_ok=True)
        if self.era == "all_eras":
            use_eras = self.all_eras
        else:
            use_eras = [self.era]
        files = [
            "/".join(
                [
                    data_dir,
                    self.era,
                    file_template.format(fold=self.branch_data["fold"]),
                ]
            )
            for file_template in self.file_templates
        ]
        # Get list of all training config targets
        allbranch_trainingconfigs = flatten_collections(
            self.input()["CreateTrainingConfig_{}".format(self.branch_data["channel"])][
                "collection"
            ]
        )
        # Filter training config targets by name in each branch
        filefilter_training_config_string = "/".join(
            [
                ".*",
                self.dir_template.format(
                    era=self.era,
                    channel=self.branch_data["channel"],
                    mass=self.branch_data["mass"],
                    batch=self.branch_data["batch_num"],
                ),
                self.file_template_config,
            ]
        )
        filefilter_training_config = re.compile(filefilter_training_config_string)
        trainingconfig = [
            target
            for target in allbranch_trainingconfigs
            if filefilter_training_config.match(target.path)
        ][0]
        # Get list of all shard targets
        allbranch_shards = flatten_collections(
            self.input()[
                "CreateTrainingDataShard_{}".format(self.branch_data["channel"])
            ]["collection"]
        )
        # Filter shard targets by name in each branch and for each process
        processes, process_classes = zip(
            *self.get_process_tuple(
                self.branch_data["channel"],
                self.branch_data["mass"],
                self.branch_data["batch_num"],
            )
        )
        shards_eras = {}
        for era in use_eras:
            filefilter_shard_strings = [
                "/".join(
                    [
                        ".*",
                        self.dir_template_in.format(
                            era=era, channel=self.branch_data["channel"]
                        ),
                        self.file_template_shard.format(process=process, fold=fold),
                    ]
                )
                for process in processes
            ]
            filefilters_shard = [
                re.compile(string) for string in filefilter_shard_strings
            ]
            shards = [
                [target for target in allbranch_shards if filefilter.match(target.path)]
                for filefilter in filefilters_shard
            ]
            shards_eras[era] = shards
        # Copy filtered shard and training config files into data directory for each era
        for era in use_eras:
            shards = [shard for shards in shards_eras[era] for shard in shards]
            for shard in shards:
                shard_name = os.path.basename(shard.path)
                shard.copy_to_local("/".join([data_dir, era, shard_name]))
        trainingconfig_name = os.path.basename(trainingconfig.path)
        trainingconfig.copy_to_local("/".join([data_dir, trainingconfig_name]))

        # self.run_command(
        #     command=["ls", "-R"],
        #     run_location=data_dir
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
                "--data-dir {}".format(data_dir),
                "--era {}".format(self.era),
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
        for file_remote, file_local in zip(self.output(), files):
            file_remote.parent.touch()
            file_remote.copy_from_local(file_local)
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
        return RunTraining(**requirements)

    def run(self):
        self.publish_message("All trainings are done!")
        self.mark_complete()
