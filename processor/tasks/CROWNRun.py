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
    sampletype = luigi.Parameter()
    era = luigi.Parameter()
    analysis = luigi.Parameter()

    def workflow_requires(self):
        requirements = super(CROWNRun, self).workflow_requires()
        requirements["datasetinfo"] = ConfigureDatasets.req(self)
        return requirements

    def requires(self):
        return {"tarball": CROWNBuild.req(self)}

    def create_branch_map(self):
        dataset = ConfigureDatasets(nick=self.nick)
        dataset.run()
        inputdata = self.input()["datasetinfo"].load()
        return {i: info for i, info in enumerate(inputdata["filelist"])}

    def output(self):
        return self.local_target("{}/ntuple_{}.root".format(self.nick, self.branch))


    def run(self):
        output = self.output()
        output.parent.touch()
        info = self.branch_data
        _inputfile = info
        _outputfile = str(output.path)
        _executable = "./{}_{}_{}".format(self.analysis, self.sampletype, self.era)
        _tarballpath = str(self.input()["tarball"].path)
        # first unpack the tarball if the exec is not there yet
        if not os.path.exists(_executable):
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
        # run CROWN
        _crown_args = [_inputfile, _outputfile]
        print("Executable: {}".format(" ".join([_executable] + _crown_args)))

        code, out, error = interruptable_popen(
            [_executable] + _crown_args, stdout=PIPE, stderr=PIPE, env=my_env
        )

        if code != 0:
            print("Error when running crown {}".format(error))
            print("Output: {}".format(out))
            print("crown returned non-zero exit status {}".format(code))
            raise Exception("crown failed")
        else:
            print("Successful")

        print("=======================================================")
