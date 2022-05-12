import luigi
import os
import subprocess
from law.util import interruptable_popen
from framework import Task
from framework import console


class CROWNBuild(Task):
    """
    Gather and compile CROWN with the given configuration
    """

    # configuration variables
    scopes = luigi.ListParameter()
    shifts = luigi.Parameter()
    build_dir = luigi.Parameter()
    install_dir = luigi.Parameter()
    eras = luigi.ListParameter()
    sampletypes = luigi.ListParameter()
    analysis = luigi.Parameter()
    production_tag = luigi.Parameter()

    env_script = os.path.join(
        os.path.dirname(__file__), "../../", "setup", "setup_crown_cmake.sh"
    )

    def output(self):
        target = self.remote_target("crown_{}_{}.tar.gz".format(self.analysis, self.production_tag)
        )
        return target

    def run(self):
        # get output file path
        output = self.output()
        # convert list to comma separated strings
        if len(self.sampletypes) == 1:
            _sampletypes = self.sampletypes[0]
        else:
            _sampletypes = ",".join(self.sampletypes)
        if len(self.eras) == 1:
            _eras = self.eras[0]
        else:
            _eras = ",".join(self.eras)
        _scopes = ",".join(self.scopes)
        _analysis = str(self.analysis)
        _shifts = str(self.shifts)
        _tag = "{}/CROWN_{}".format(self.production_tag, _analysis)
        _install_dir = os.path.join(str(self.install_dir), _tag)
        _build_dir = os.path.join(str(self.build_dir), _tag)
        _crown_path = os.path.abspath("CROWN")
        _compile_script = os.path.join(
            str(os.path.abspath("processor")), "tasks", "compile_crown.sh"
        )

        if os.path.exists(output.path):
            console.log("tarball already existing in {}".format(output.path))

        elif os.path.exists(os.path.join(_install_dir, output.basename)):
            console.log(
                "tarball already existing in tarball directory {}".format(_install_dir)
            )
            output.copy_from_local(os.path.join(_install_dir, output.basename))
        else:
            console.rule("Building new CROWN tarball")
            # create build directory
            if not os.path.exists(_build_dir):
                os.makedirs(_build_dir)
            _build_dir = os.path.abspath(_build_dir)
            # same for the install directory
            if not os.path.exists(_install_dir):
                os.makedirs(_install_dir)
            _install_dir = os.path.abspath(_install_dir)

            # set environment variables
            my_env = self.set_environment(self.env_script)

            # checking cmake path
            code, _cmake_executable, error = interruptable_popen(
                ["which", "cmake"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=my_env,
            )
            # actual payload:
            console.rule("Starting cmake step for CROWN")
            console.log("Using cmake {}".format(_cmake_executable.replace("\n", "")))
            console.log("Using CROWN {}".format(_crown_path))
            console.log("Using build_directory {}".format(_build_dir))
            console.log("Using install directory {}".format(_install_dir))
            console.log("Settings used: ")
            console.log("Analysis: {}".format(_analysis))
            console.log("Sampletype: {}".format(_sampletypes))
            console.log("Era: {}".format(_eras))
            console.log("Channels: {}".format(_scopes))
            console.log("Shifts: {}".format(_shifts))
            console.rule("")

            # run crown compilation script
            command = [
                "bash",
                _compile_script,
                _crown_path,  # CROWNFOLDER=$1
                _analysis,  # ANALYSIS=$2
                _sampletypes,  # SAMPLES=$3
                _eras,  # ERAS=$4
                _scopes,  # SCOPES=$5
                _shifts,  # SHIFTS=$6
                _install_dir,  # INSTALLDIR=$7
                _build_dir,  # BUILDDIR=$8
                output.basename,  # TARBALLNAME=$9
            ]
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                bufsize=1,
                universal_newlines=True,
            ) as p:
                for line in p.stdout:
                    if line != "\n":
                        console.log(line.replace("\n", ""))

            if p.returncode != 0:
                console.log("Error when building crown {}".format(command))
                console.log("Output: {}".format(p.stderr))
                console.log(
                    "building returned non-zero exit status {}".format(p.returncode)
                )
                console.rule()
                raise Exception("crown build failed")
            else:
                console.rule("Successful crown build ! ")
            console.log(
                "Copying from local: {}".format(
                    os.path.join(_install_dir, output.basename)
                )
            )
            output.parent.touch()
            console.log("Copying to remote: {}".format(output.path))
            output.copy_from_local(os.path.join(_install_dir, output.basename))
        console.rule("Finished CROWNRun")
