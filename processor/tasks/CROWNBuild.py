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

        if os.path.exists(output.path):
            console.log("tarball already existing in {}".format(output.path))

        # elif os.path.exists(os.path.join(_install_dir, output.basename)):
        #     console.log(
        #         "tarball already existing in tarball directory {}".format(_install_dir)
        #     )
        #     output.copy_from_local(os.path.join(_install_dir, output.basename))
        else:
            console.log("Building new tarball")
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

            # run CROWN build step
            _cmake_cmd = [_cmake_executable.replace("\n", ""), _crown_path]

            _cmake_args = [
                "-DANALYSIS={ANALYSIS}".format(ANALYSIS=_analysis),
                "-DSAMPLES={SAMPLES}".format(SAMPLES=_sampletype),
                "-DERAS={ERAS}".format(ERAS=_era),
                "-DCHANNELS={CHANNELS}".format(CHANNELS=_channels),
                "-DSHIFTS={SHIFTS}".format(SHIFTS=_shifts),
                "-DINSTALLDIR={INSTALLDIR}".format(INSTALLDIR=_install_dir),
                "-B{BUILDFOLDER}".format(BUILDFOLDER=_build_dir),
            ]
            console.rule("Running cmake")
            console.log("{}".format(" ".join(_cmake_cmd + _cmake_args)))
            console.rule()

            code, out, error = interruptable_popen(
                _cmake_cmd + _cmake_args, stdout=PIPE, stderr=PIPE, env=my_env
            )
            for stdout_line in iter(code.stdout.readline, ""):
                yield stdout_line
            # if successful save Herwig-cache and run-file as tar.gz
            if code != 0:
                console.log("Error when running cmake {}".format(error))
                console.log("Output: {}".format(out))
                console.log("cmake returned non-zero exit status {}".format(code))
                console.rule()
                raise Exception("cmake failed")
            else:
                console.rule("Successful cmake build !")

            console.rule("Running make")
            console.log("{}".format(" ".join(["make", "install"])))
            console.rule()
            code, out, error = interruptable_popen(
                ["make", "install"],
                stdout=PIPE,
                stderr=PIPE,
                env=my_env,
                cwd=_build_dir,
            )
            if code != 0:
                console.log("Error when running make {}".format(error))
                console.log("Output: {}".format(out))
                console.log("make returned non-zero exit status {}".format(code))
                console.rule()
                raise Exception("make failed")
            else:
                console.rule("Successful build !")

            code, out, error = interruptable_popen(
                ["touch", output.basename],
                stdout=PIPE,
                stderr=PIPE,
                env=my_env,
                cwd=os.path.join(_install_dir),
            )
            command = [
                "tar",
                "-czvf",
                output.basename,
                "--exclude={}".format(output.basename),
                ".",
            ]
            console.rule("Running tar")
            console.log("{}".format(" ".join(command)))
            console.rule()
            code, out, error = interruptable_popen(
                command,
                stdout=PIPE,
                stderr=PIPE,
                env=my_env,
                cwd=os.path.join(_install_dir),
            )
            if code != 0:
                console.log("Error when creating tarball {}".format(error))
                console.log("Output: {}".format(out))
                console.log("tar returned non-zero exit status {}".format(code))
                console.rule()
                raise Exception("tar failed")
            else:
                console.rule("Successful tarball creation ! ")
            output.copy_from_local(os.path.join(_install_dir, output.basename))
