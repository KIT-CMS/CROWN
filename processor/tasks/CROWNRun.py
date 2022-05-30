import law
import luigi
import os
from CROWNBuild import CROWNBuild
import tarfile
from ConfigureDatasets import ConfigureDatasets
import subprocess
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

    output_collection_cls = law.NestedSiblingFileCollection

    nick = luigi.Parameter()
    scopes = luigi.ListParameter()
    sampletype = luigi.Parameter()
    sampletypes = luigi.ListParameter()
    era = luigi.Parameter()
    eras = luigi.ListParameter()
    analysis = luigi.Parameter()
    config = luigi.Parameter()
    production_tag = luigi.Parameter()
    files_per_task = luigi.IntParameter(default=1)

    def htcondor_job_config(self, config, job_num, branches):
        config = super().htcondor_job_config(config, job_num, branches)
        config.custom_content.append(
            ("JobBatchName", f"{self.nick}-{self.analysis}-{self.config}")
        )
        return config

    def modify_polling_status_line(self, status_line):
        """
        Hook to modify the status line that is printed during polling.
        """
        return f"{status_line} - {law.util.colored(self.nick, color='light_cyan')}"

    def workflow_requires(self):
        requirements = super(CROWNRun, self).workflow_requires()
        requirements["datasetinfo"] = ConfigureDatasets.req(self)
        requirements["tarball"] = CROWNBuild.req(self)
        return requirements

    def requires(self):
        return {"tarball": CROWNBuild.req(self)}

    def create_branch_map(self):
        dataset = ConfigureDatasets(nick=self.nick, production_tag=self.production_tag)
        dataset.run()
        with self.input()["datasetinfo"].localize("r") as _file:
            inputdata = _file.load()
        branch_map = {}
        if len(inputdata["filelist"]) == 0:
            raise Exception("No files found for dataset {}".format(self.nick))
        for i, filename in enumerate(inputdata["filelist"]):
            if (int(i / self.files_per_task)) not in branch_map:
                branch_map[int(i / self.files_per_task)] = []
            branch_map[int(i / self.files_per_task)].append(filename)
        return branch_map

    def output(self):
        nicks = [
            "{era}/{nick}/{scope}/{nick}_{branch}.root".format(
                era=self.era,
                nick=self.nick,
                branch=self.branch,
                scope=scope,
            )
            for scope in self.scopes
        ]
        targets = self.remote_targets(nicks)
        for target in targets:
            target.parent.touch()
        return targets

    def run(self):
        outputs = self.output()
        info = self.branch_data
        _workdir = os.path.abspath("workdir")
        ensure_dir(_workdir)
        _inputfiles = info
        # set the outputfilename to the first name in the output list, removing the scope suffix
        _outputfile = str(
            outputs[0].basename.replace("_{}.root".format(self.scopes[0]), ".root")
        )
        _executable = "{}/{}_{}_{}".format(
            _workdir, self.config, self.sampletype, self.era
        )
        console.log(
            "Getting CROWN tarball from {}".format(self.input()["tarball"].uri())
        )
        with self.input()["tarball"].localize("r") as _file:
            _tarballpath = _file.path
        # first unpack the tarball if the exec is not there yet
        if os.path.exists(
            "unpacking_{}_{}_{}".format(self.config, self.sampletype, self.era)
        ):
            time.sleep(5)
        if not os.path.exists(_executable):
            open(
                "unpacking_{}_{}_{}".format(self.config, self.sampletype, self.era),
                "a",
            ).close()
            tar = tarfile.open(_tarballpath, "r:gz")
            tar.extractall("workdir")
            os.remove(
                "unpacking_{}_{}_{}".format(self.config, self.sampletype, self.era)
            )
        # set environment using env script
        my_env = self.set_environment("{}/init.sh".format(_workdir))
        _crown_args = [_outputfile] + _inputfiles
        _executable = "./{}_{}_{}".format(self.config, self.sampletype, self.era)
        # actual payload:
        console.rule("Starting CROWNRun")
        console.log("Executable: {}".format(_executable))
        console.log("inputfile {}".format(_inputfiles))
        console.log("outputfile {}".format(_outputfile))
        console.log("workdir {}".format(_workdir))  # run CROWN
        with subprocess.Popen(
            [_executable] + _crown_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            bufsize=1,
            universal_newlines=True,
            env=my_env,
            cwd=_workdir,
        ) as p:
            for line in p.stdout:
                if line != "\n":
                    console.log(line.replace("\n", ""))
            for line in p.stderr:
                if line != "\n":
                    console.log("Error: {}".format(line.replace("\n", "")))
        if p.returncode != 0:
            console.log(
                "Error when running crown {}".format(
                    [_executable] + _crown_args,
                )
            )
            console.log("crown returned non-zero exit status {}".format(p.returncode))
            raise Exception("crown failed")
        else:
            console.log("Successful")
        console.log("Output files afterwards: {}".format(os.listdir(_workdir)))
        for i, outputfile in enumerate(outputs):
            outputfile.parent.touch()
            # for each outputfile, add the scope suffix
            outputfile.copy_from_local(
                os.path.join(
                    _workdir,
                    _outputfile.replace(".root", "_{}.root".format(self.scopes[i])),
                )
            )
        console.rule("Finished CROWNRun")
