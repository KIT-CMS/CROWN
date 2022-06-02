# coding: utf-8

# """
# Simple law tasks that demonstrate how use remote file targets and how to run on remote resources using HTCondor.
# """


import sys
import os
import luigi
import law
from shutil import rmtree
from framework import Task, HTCondorWorkflow
from law import LocalWorkflow
from tempfile import mkdtemp
from CreateTrainingDatasets import CreateTrainingDatasets, CreateTrainingDatasetsAllEras
from copy import copy
law.contrib.load("tasks")  # to have the RunOnceTask

# law.contrib.load("tasks")  # to have the RunOnceTask

#Task runs over HTCondor
class RunTraining(HTCondorWorkflow, law.LocalWorkflow):
    era = luigi.Parameter(description="Run era")
    channel = luigi.Parameter(description="Decay Channel")
    mass = luigi.Parameter(description="Mass hypothesis of heavy NMSSM Higgs boson.")
    batch_num = luigi.Parameter(description="Group of mass hypotheses of light NMSSM Higgs boson.")

    files_template = [
        "{prefix}fold{fold}_keras_model.h5",
        "{prefix}fold{fold}_keras_preprocessing.pickle",
        "{prefix}fold{fold}_loss.pdf",
        "{prefix}fold{fold}_loss.png",
    ]
        # "{prefix}training{fold}.log"

    def create_branch_map(self):
        return [
            {
                "era": self.era, 
                "channel": self.channel, 
                "mass": self.mass, 
                "batch": self.batch_num, 
                "fold":fold
            } 
            for fold in ["0","1"]
        ]

    # Requirements dependant on whether all_eras was used
    def requires(self):
        return self.workflow_requires()

            
    def workflow_requires(self):
        requirements = super(RunTraining, self).workflow_requires()
        requirements_args = {
            "era": self.era, 
            "channel": self.channel, 
            "mass": self.mass, 
            "batch_num": self.batch_num
        }
        # if "default/" in self.production_tag:
        #     requirements_args["production_tag"] = self.production_tag
        if self.era == "all_eras":
            requirements["CreateTrainingDatasets"] = CreateTrainingDatasets(**requirements_args)
            requirements["CreateTrainingDatasetsAllEras"] = CreateTrainingDatasetsAllEras(**requirements_args)
        else:
            requirements["CreateTrainingDatasets"] = CreateTrainingDatasets(**requirements_args)
        return requirements



    def output(self):
        # Define output files: fold 0 and 1 of training data + config file
        prefix = "training_dataset_{era}_{channel}_{mass}_{batch}/".format(
                    era=self.era,
                    channel=self.channel,
                    mass=self.mass,
                    batch=self.batch_num
                )
        
        fold = self.branch_data["fold"]
        files = [file_string.format(prefix=prefix, fold=fold) for file_string in self.files_template]
        targets = self.remote_targets(files)
        for target in targets:
            target.parent.touch()
        return targets

    def run(self):
        cmb_tag = "{channel}_{mass}_{batch}".format(
                    channel=self.channel,
                    mass=self.mass,
                    batch=self.batch_num
                )
        fold = self.branch_data["fold"]
        # self.publish_message("This is layer 0: {}".format(self.input()))
        # self.publish_message("This is layer 1a: {}".format(self.input()['CreateTrainingDatasets']['collection']))
        # self.publish_message("This is layer 2a: {}".format(self.input()['CreateTrainingDatasets']['collection'][0]))
        # self.publish_message("This is layer 1b: {}".format(self.input()['CreateTrainingDatasetsAllEras']))
        collections = self.input()['CreateTrainingDatasets']['collection']
        prefix = "sm-htt-analysis/output/ml/{era}_{cmb_tag}/".format(era=self.era, cmb_tag=cmb_tag)
        if self.era == "all_eras":
            all_eras_dataset_config = self.input()['CreateTrainingDatasetsAllEras']
            eras = ["2016", "2017", "2018"]
            for i_era, era in enumerate(eras):
                prefix_loop = "sm-htt-analysis/output/ml/{era}_{cmb_tag}/".format(era=era, cmb_tag=cmb_tag)
                for i_file, copyFile in enumerate(collections[i_era]):
                    if i_file==int(fold):
                        print("copying {}".format(copyFile))
                        filename = os.path.basename(copyFile.path)
                        copyFile.copy_to_local(prefix_loop + filename)
            print("copying {}".format(all_eras_dataset_config))
            filename = os.path.basename(all_eras_dataset_config.path)
            all_eras_dataset_config.copy_to_local(prefix + filename)
        else:
            for i_file, copyFile in enumerate(collections[0]):
                if (i_file==2) or (i_file==int(fold)):
                    print("copying {}".format(copyFile))
                    filename = os.path.basename(copyFile.path)
                    copyFile.copy_to_local(prefix + filename)

        # self.run_command(
        #     command=["ls", "-R"],
        #     run_location=prefix
        # )
        self.run_command(
            command=[
                "ml/run_training.sh", 
                self.era, 
                self.channel, 
                "{}_{}".format(self.mass, self.batch_num),
                fold
            ], 
            run_location="sm-htt-analysis"
        )
        # self.run_command(
        #     command=["ls", "-R"],
        #     run_location=prefix
        # )
        files = [file_string.format(prefix=prefix, fold=fold) for file_string in self.files_template]
        for file_remote, file_local in zip(self.output(), files):
            file_remote.copy_from_local(file_local)