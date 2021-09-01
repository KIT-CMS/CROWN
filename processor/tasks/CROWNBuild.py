from platform import processor
import law
import luigi
import os
from subprocess import PIPE
from law.util import interruptable_popen
from framework import Task
from framework import console


class CROWNBuild(Task):
    """
    Gather and compile CROWN with the given configuration
    """

    # configuration variables
    channels = luigi.Parameter()
    shifts = luigi.Parameter()
    build_dir = luigi.Parameter()
    install_dir = luigi.Parameter()
    era = luigi.Parameter()
    sampletype = luigi.Parameter()
    analysis = luigi.Parameter()
    production_tag = luigi.Parameter()

    env_script = os.path.join(
        os.path.dirname(__file__), "../../", "setup", "setup_crown_cmake.sh"
    )

    def output(self):
        target = self.remote_target(
            "{}/crown_{}_{}.tar.gz".format(
                self.production_tag, self.era, self.sampletype
            )
        )
        target.parent.touch()
        return target

    def run(self):
        # get output file path
        output = self.output()
        output.parent.touch()
        _sampletype = str(self.sampletype)
        _era = str(self.era)
        _channels = str(self.channels)
        _analysis = str(self.analysis)
        _shifts = str(self.shifts)
        _tag = "{}_{}".format(_era, _sampletype)
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
            console.rule("Building new tarball")
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
                ["which", "cmake"], stdout=PIPE, stderr=PIPE, env=my_env
            )
            # actual payload:
            console.rule("Starting cmake step for CROWN")
            console.log("Using cmake {}".format(_cmake_executable.replace("\n", "")))
            console.log("Using CROWN {}".format(_crown_path))
            console.log("Using build_directory {}".format(_build_dir))
            console.log("Using install directory {}".format(_install_dir))
            console.log("Settings used: ")
            console.log("Analysis: {}".format(_analysis))
            console.log("Sampletype: {}".format(_sampletype))
            console.log("Era: {}".format(_era))
            console.log("Channels: {}".format(_channels))
            console.log("Shifts: {}".format(_shifts))
            console.rule("")

            # run crown compilation script
            command = [
                "bash",
                _compile_script,
                _crown_path,  # CROWNFOLDER=$1
                _analysis,  # ANALYSIS=$2
                _sampletype,  # SAMPLES=$3
                _era,  # ERAS=$4
                _channels,  # CHANNEL=$5
                _shifts,  # SHIFTS=$6
                _install_dir,  # INSTALLDIR=$7
                _build_dir,  # BUILDDIR=$8
                output.basename,  # TARBALLNAME=$9
            ]
            code, out, error = interruptable_popen(
                command,
                rich_console=console,
                stdout=PIPE,
                stderr=PIPE,
            )
            if code != 0:
                console.log("Error when building crown {}".format(error))
                console.log("Output: {}".format(out))
                console.log("tar returned non-zero exit status {}".format(code))
                console.rule()
                raise Exception("crown build failed")
            else:
                console.rule("Successful crown build ! ")
            output.copy_from_local(os.path.join(_install_dir, output.basename))
