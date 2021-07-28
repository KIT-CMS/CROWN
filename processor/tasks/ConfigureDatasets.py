import law
import luigi
import os

from subprocess import PIPE
from law.util import interruptable_popen

from processor.framework import Task


class ConfigureDatasets(Task):
    """
    Gather information on the selected datasets, for now just mentioning an explicit dataset
    """
    datasets = luigi.Parameter()

    def output(self):
        return self.local_target("")


    def run(self):
        _datasets = str(self.datasets)
