# coding: utf-8

# """
# Simple law tasks that demonstrate how use remote file targets and how to run on remote resources using HTCondor.
# """


import os
import luigi
import law
from shutil import rmtree
from framework import Task, HTCondorWorkflow
from tempfile import mkdtemp

# law.contrib.load("tasks")  # to have the RunOnceTask

#Task runs over HTCondor
class CreateTrainingDatasets(HTCondorWorkflow, law.LocalWorkflow):
    
    era = luigi.Parameter(description="Run era")
    channel = luigi.Parameter(description="Decay Channel")
    mass = luigi.Parameter(description="Mass hypothesis of heavy NMSSM Higgs boson.")
    batch_num = luigi.Parameter(description="Group of mass hypotheses of light NMSSM Higgs boson.")

    files_template = [
        "{prefix}fold0_training_dataset.root",
        "{prefix}fold1_training_dataset.root",
        "{prefix}dataset_config.yaml",
    ]
    def create_branch_map(self):
        if self.era == "all_eras":
            eras = ["2016", "2017", "2018"]
        else:
            eras = [self.era]
        return [{"era": era, "channel": self.channel, "mass": self.mass, "batch": self.batch_num} for era in eras]

    def output(self):
        # Define output files: fold 0 and 1 of training data + config file
        prefix = "training_dataset_{era}_{channel}_{mass}_{batch}/".format(
                    era=self.branch_data["era"],
                    channel=self.channel,
                    mass=self.mass,
                    batch=self.batch_num
                )
        files = [file_string.format(prefix=prefix) for file_string in self.files_template]
        targets = self.remote_targets(files)
        for target in targets:
            target.parent.touch()
        return targets

    def run(self):
        # Create training datasets
        self.run_command(
            command=[
                "ml/create_training_dataset.sh", 
                self.branch_data["era"], 
                self.channel, 
                self.mass, 
                self.batch_num
            ],
            run_location="sm-htt-analysis"
        )
        # Copy resulting files to remote storage
        prefix = "sm-htt-analysis/output/ml/{era}_{channel}_{mass}_{batch}/".format(
                    era=self.branch_data["era"],
                    channel=self.channel,
                    mass=self.mass,
                    batch=self.batch_num
                )
        self.run_command(
            command=["ls", "-R"],
            run_location=prefix
        )
        files = [file_string.format(prefix=prefix) for file_string in self.files_template]
        for file_remote, file_local in zip(self.output(), files):
            file_remote.copy_from_local(file_local)

# Task runs local
class CreateTrainingDatasetsAllEras(Task):
    era = luigi.Parameter(description="Run era")
    channel = luigi.Parameter(description="Decay Channel")
    mass = luigi.Parameter(description="Mass hypothesis of heavy NMSSM Higgs boson.")
    batch_num = luigi.Parameter(description="Group of mass hypotheses of light NMSSM Higgs boson.")

    # Requirements dependant on whether all_eras was used
    def requires(self):
        if self.era!="all_eras":
            raise Exception(
                "CreateTrainingDatasetsAllEras task is intended for all_eras, but {} was given.".format(
                    self.era
                )
            )
        requirements_args = {
            "era": self.era, 
            "channel": self.channel, 
            "mass": self.mass, 
            "batch_num": self.batch_num
        }
        return CreateTrainingDatasets(**requirements_args)


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
            # Fetch config files of all eras to local
            for i, era in enumerate(["2016", "2017", "2018"]):
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
            self.run_command(command=command, run_location="sm-htt-analysis")
            prefix = "{tmpdir}/all_eras_{cmb_tag}/".format(tmpdir=temporary_dir, cmb_tag=cmb_tag)
            # Send combined configs to remote
            self.output().touch()
            self.output().copy_from_local("{prefix}dataset_config.yaml".format(prefix=prefix))