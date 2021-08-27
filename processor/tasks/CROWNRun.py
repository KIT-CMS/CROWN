import law
import luigi
import os
from CROWNBuild import CROWNBuild
import tarfile
from ConfigureDatasets import ConfigureDatasets
from subprocess import PIPE
from law.util import interruptable_popen
import time
from framework import console

from framework import Task, HTCondorWorkflow


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


class CROWNRun(Task, HTCondorWorkflow, law.LocalWorkflow):
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
        requirements["tarball"] = CROWNBuild.req(self)
        return requirements

    def requires(self):
        return {"tarball": CROWNBuild.req(self)}

    def create_branch_map(self):
        dataset = ConfigureDatasets(nick=self.nick)
        dataset.run()
        with self.input()["datasetinfo"].localize('r') as _file:
            inputdata = _file.load()
        return {i: info for i, info in enumerate(inputdata["filelist"])}

    def output(self):
        return self.remote_target("{}/ntuple_{}.root".format(self.nick, self.branch))

    def run(self):
        output = self.output()
        output.parent.touch()
        info = self.branch_data
        _workdir = os.path.abspath("workdir")
        ensure_dir(_workdir)
        _inputfile = info
        _outputfile = str(output.path).replace(".root", "_running.root")
        _executable = "{}/{}_{}_{}".format(
            _workdir, self.analysis, self.sampletype, self.era
        )
        with self.input()['tarball'].localize('r') as _file:
            _tarballpath = _file.path
        # first unpack the tarball if the exec is not there yet
        if os.path.exists(
            "unpacking_{}_{}_{}".format(self.analysis, self.sampletype, self.era)
        ):
            time.sleep(5)
        if not os.path.exists(_executable):
            open(
                "unpacking_{}_{}_{}".format(self.analysis, self.sampletype, self.era),
                "a",
            ).close()
            tar = tarfile.open(_tarballpath, "r:gz")
            tar.extractall("workdir")
            os.remove(
                "unpacking_{}_{}_{}".format(self.analysis, self.sampletype, self.era)
            )
        # set environment using env script
        my_env = self.set_environment("{}/init.sh".format(_workdir))
        _crown_args = [_inputfile, _outputfile]
        _executable = "./{}_{}_{}".format(self.analysis, self.sampletype, self.era)
        # actual payload:
        console.rule("Starting CROWNRun")
        console.log("Executable: {}".format(_executable))
        console.log("inputfile {}".format(_inputfile))
        console.log("outputfile {}".format(_outputfile))
        console.log("workdir {}".format(_workdir))  # run CROWN
        code, out, error = interruptable_popen(
            [_executable] + _crown_args,
            stdout=PIPE,
            stderr=PIPE,
            env=my_env,
            cwd=_workdir,
        )
        if code != 0:
            console.log("Error when running crown {}".format(error))
            console.log("Output: {}".format(out))
            console.log("crown returned non-zero exit status {}".format(code))
            raise Exception("crown failed")
        else:
            console.log("Successful")
            console.log("Output: {}".format(out))
        output.move_from_local(_outputfile)
        console.rule("Finished CROWNRun")
