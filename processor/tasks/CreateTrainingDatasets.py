# coding: utf-8

# """
# Simple law tasks that demonstrate how use remote file targets and how to run on remote resources using HTCondor.
# """


import os
import luigi
import law
from shutil import rmtree
from framework import Task, HTCondorWorkflow
from law import LocalWorkflow
from tempfile import mkdtemp

# law.contrib.load("tasks")  # to have the RunOnceTask

#Task runs over HTCondor
class CreateTrainingDatasets(HTCondorWorkflow, LocalWorkflow):

# class CreateTrainingDatasetsShard(Task, LocalWorkflow):
    era = luigi.Parameter()
    channel = luigi.Parameter()
    mass = luigi.Parameter()
    batch_num = luigi.Parameter()

    def create_branch_map(self):
        if self.era == "all_eras":
            eras = ["2016", "2017", "2018"]
        else:
            eras = [self.era]
        return {i: data for i, data in enumerate(
            [{"era": era, "channel": self.channel, "mass": self.mass, "batch": self.batch_num} for era in eras]
            )}

    def output(self):
        # Define output files: fold 0 and 1 of training data + config file
        prefix = "training_dataset_{era}_{channel}_{mass}_{batch}/".format(
                    era=self.branch_data["era"],
                    channel=self.channel,
                    mass=self.mass,
                    batch=self.batch_num
                )
        files = [
            "{prefix}fold0_training_dataset.root".format(prefix=prefix),
            "{prefix}fold1_training_dataset.root".format(prefix=prefix),
            "{prefix}dataset_config.yaml".format(prefix=prefix)
        ]
        targets = self.remote_targets(files)
        targets[0].parent.touch()
        return targets

    def run(self):
        # Create training datasets
        self.run_command(command=["ml/create_training_dataset.sh", self.branch_data["era"], self.channel, self.mass, self.batch_num], run_location="ml_scripts")
        # Copy resulting files to remote storage
        prefix = "ml_scripts/output/ml/{era}_{channel}_{mass}_{batch}/".format(
                    era=self.branch_data["era"],
                    channel=self.channel,
                    mass=self.mass,
                    batch=self.batch_num
                )
        self.output()[0].copy_from_local("{prefix}fold0_training_dataset.root".format(prefix=prefix))
        self.output()[1].copy_from_local("{prefix}fold1_training_dataset.root".format(prefix=prefix))
        self.output()[2].copy_from_local("{prefix}dataset_config.yaml".format(prefix=prefix))

# Task runs local
class CreateTrainingDatasetsAllEras(Task):
    era = luigi.Parameter()
    channel = luigi.Parameter()
    mass = luigi.Parameter()
    batch_num = luigi.Parameter()

    # Requirements dependant on whether all_eras was used
    def requires(self):
        if self.era!="all_eras":
            raise Exception(
                "CreateTrainingDatasetsAllEras task is intended for all_eras, but {} was given.".format(
                    self.era
                )
            )
        return CreateTrainingDatasets.req(self)


    def output(self):
        # Require combined config file
        prefix = "training_dataset_all_eras_{channel}_{mass}_{batch}/".format(
                    channel=self.channel,
                    mass=self.mass,
                    batch=self.batch_num
                )
        files = "{prefix}dataset_config.yaml".format(prefix=prefix)
        target = self.remote_target(files)
        target.parent.touch()
        return target

    def run(self):
        cmb_tag = "{channel}_{mass}_{batch}".format(
                    channel=self.channel,
                    mass=self.mass,
                    batch=self.batch_num
                )
        # If all_eras:
        if self.era == "all_eras":
            temporary_dir = mkdtemp(dir="/tmp/{user}".format(user=self.user_name))
            print(temporary_dir)
            # print(self.input()["collection"])
            # Fetch config files of all eras to local
            for i, era in enumerate(["2016", "2017", "2018"]):
                # print(self.input()["collection"][i][2])
                prefix_ = "{tmpdir}/{era}_{cmb_tag}/".format(tmpdir=temporary_dir, era=era, cmb_tag=cmb_tag)
                era_conf_file = self.input()["collection"][i][2] #dataset_config.yaml
                era_conf_file.copy_to_local("{prefix}dataset_config.yaml".format(prefix=prefix_))
            # Combine configs
            command = ["ml/combine_configs.sh", 
                self.channel, 
                self.mass, 
                self.batch_num,
                temporary_dir
            ]
            self.run_command(command=command, run_location="ml_scripts")
            prefix = "{tmpdir}/all_eras_{cmb_tag}/".format(tmpdir=temporary_dir, cmb_tag=cmb_tag)
            # Send combined configs to remote
            self.output().touch()
            self.output().copy_from_local("{prefix}dataset_config.yaml".format(prefix=prefix))