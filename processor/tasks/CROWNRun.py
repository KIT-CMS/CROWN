import law
import luigi
import os
from CROWNBuild import CROWNBuild
import tarfile
from ConfigureDatasets import ConfigureDatasets
from subprocess import PIPE
from law.util import interruptable_popen

from framework import Task, HTCondorWorkflow


class CROWNRun(Task, law.LocalWorkflow):
    """
    Gather and compile CROWN with the given configuration
    """

    output_collection_cls = law.SiblingFileCollection

    nick = luigi.Parameter()

    def workflow_requires(self):
        return {"datasetinfo": ConfigureDatasets.req(self, self.nick)}

    def create_branch_map(self):
        print(self.input())
        datasets = self.input()["datasetinfo"].load()
        return {i: info for i, info in enumerate(datasets["filelist"])}

    def output(self):
        return self.local_target("ntuple_{}.root".format(self.branch))

    def requires(self):
        return {"tarball": CROWNBuild.req(self)}

    def run(self):
        output = self.output()
        info = self.branch_data
        _inputfile = info
        _outputfile = str(output.path)
        _tarballpath = str(self.input()["tarball"].path)
        # first unpack the tarball
        tar = tarfile.open(_tarballpath, "r:gz")
        tar.extractall()

        # set environment using env script
        my_env = self.set_environment("./init.sh")

        # actual payload:
        print("=========================================================")
        print("| Starting CROWN")
        print("| inputfile {}".format(_inputfile))
        print("| outputfile {}".format(_outputfile))
        print("=========================================================")
        # run CROWN build step
        _crown_cmd = ["./{}_{}_{}".format(self.analysis, self.sampletype, self.era)]

        _crown_args = [_inputfile, _outputfile]
        print("Executable: {}".format(" ".join(_crown_cmd + _crown_args)))

        code, out, error = interruptable_popen(
            _crown_cmd + _crown_args, stdout=PIPE, stderr=PIPE, env=my_env
        )

        if code != 0:
            print("Error when running crown {}".format(error))
            print("Output: {}".format(out))
            print("crown returned non-zero exit status {}".format(code))
            raise Exception("crown failed")
        else:
            print("Successful")

        print("=======================================================")
