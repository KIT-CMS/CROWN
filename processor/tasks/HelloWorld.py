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
class LawBase(law.Task):
    def local_path(self, *path):
        # LAW_DATA_PATH is defined in setup.sh
        parts = (os.getenv("ANALYSIS_DATA_PATH"),) + path
        return os.path.join(*(str(p) for p in parts))

    def local_target(self, *path, **kwargs):
        return law.LocalFileTarget(self.local_path(*path), **kwargs)


class HelloWorldSave(Task):
    def output(self):
        target = self.remote_target("text.txt")
        target.parent.touch()
        return target

    def run(self):
        saveText = "This is "
        self.output().dump(saveText)

class HelloWorldPrint(Task, HTCondorWorkflow, law.LocalWorkflow):
    def requires(self):
        return HelloWorldSave.req(self)

    def output(self):
        target = self.remote_target("whole.txt")
        target.parent.touch()
        return target

    def create_branch_map(self):
        return [0]

    def run(self):
        import tensorflow as tf
        tf.test.is_gpu_available()
        readText = self.input().load()
        self.publish_message("This is the input: {}".format(self.input()))
        wholeText = readText + "a triumph!"

        self.publish_message("This is the text: {}".format(wholeText))
        output = self.output()
        output.parent.touch()
        output.dump(wholeText)


class HelloWorldPrintFinal(Task, law.tasks.RunOnceTask):
    def requires(self):
        return HelloWorldPrint.req(self)
    
    def run(self):
        self.publish_message("This is the 0: {}".format(self.input()))
        self.publish_message("This is the 1: {}".format(self.input()['collection']))
        self.publish_message("This is the 2: {}".format(self.input()['collection'][0]))
        self.publish_message("This is the 3: {}".format(self.input()['collection'][0].load()))
        self.mark_complete()