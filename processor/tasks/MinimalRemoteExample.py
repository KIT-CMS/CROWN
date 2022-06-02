# # coding: utf-8

# """
# Simple law tasks that demonstrate how to build up a task tree with some outputs and dependencies.

# The first task (FetchLoremIpsum) will download 1 of 6 different versions of a "lorem ipsum" text.
# The next task (CountChars) determines and saves the frequency of every character in a json file.
# After that, the count files are merged (MergeCounts). The last task (ShowFrequencies) illustrates
# the "measured" frequencies and prints the result which is also sent as a message to the scheduler.
# """


import os
import luigi
import law
from framework import Task, HTCondorWorkflow

law.contrib.load("tasks")  # to have the RunOnceTask

class SaveToRemote(Task):
    # Output target is remote and accessed using gfal2.
    def output(self):
        target = self.remote_target("RemoteFileIn.txt")
        target.parent.touch()
        return target

    def run(self):
        saveText = "This is "
        self.output().dump(saveText)

class RunRemote(HTCondorWorkflow, law.LocalWorkflow):
    def requires(self):
        return SaveToRemote.req(self)

    def workflow_requires(self):
        requirements = super(RunRemote, self).workflow_requires()
        requirements["TextFile"] = SaveToRemote.req(self)
        return requirements

    # Output target is remote and accessed using gfal2.
    def output(self):
        target = self.remote_target("RemoteFileChanged.txt")
        target.parent.touch()
        return target

    # Branch map is necessary for remote tasks.
    def create_branch_map(self):
        return [0]

    # Task is run remotely and has access to GPU resources.
    def run(self):
        try:
            import tensorflow as tf
            tf.test.is_gpu_available()
        except:
            print("Tensorflow not found.")
        readText = self.input().load()
        self.publish_message("This is the input: {}".format(self.input()))
        wholeText = readText + "a triumph!"

        self.publish_message("This is the text: {}".format(wholeText))
        output = self.output()
        output.parent.touch()
        output.dump(wholeText)


class ReadFromRemote(Task, law.tasks.RunOnceTask):
    def requires(self):
        return RunRemote.req(self)
    
    # Print all layers of the input individually. Final layer is the requested data.
    def run(self):
        self.publish_message("This is layer 0: {}".format(self.input()))
        self.publish_message("This is layer 1: {}".format(self.input()['collection']))
        self.publish_message("This is layer 2: {}".format(self.input()['collection'][0]))
        self.publish_message("This is layer 3: {}".format(self.input()['collection'][0].load()))
        self.mark_complete()