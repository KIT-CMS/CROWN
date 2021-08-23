import law
import luigi
import os

from subprocess import PIPE
from law.util import interruptable_popen

from framework import Task


class CROWNBuild(Task):
    """
    Gather and compile CROWN with the given configuration
    """

    # configuration variables
    channels = luigi.Parameter()
    analysis = luigi.Parameter()
    shifts = luigi.Parameter()
    build_dir = luigi.Parameter()
    install_dir = luigi.Parameter()
    env_script = os.path.join(
        os.path.dirname(__file__), "../../", "setup", "setup_crown_cmake.sh"
    )

    def output(self):
        return self.local_target("crown_{}_{}.tar.gz".format(self.era, self.sampletype))

    def run(self):
        _sampletype = str(self.sampletype)
        _era = str(self.era)
        _channels = str(self.channels)
        _analysis = str(self.analysis)
        _shifts = str(self.shifts)
        _tag = "{}_{}".format(_era, _sampletype)

        _build_dir = os.path.join(str(self.build_dir), _tag)
        _install_dir = os.path.join(str(self.install_dir), _tag)
        # find crown
        _crown_path = os.path.abspath("CROWN")

        # create build directory
        if not os.path.exists(_build_dir):
            os.makedirs(_build_dir)
        _build_dir = os.path.abspath(_build_dir)
        # same for the install directory
        if not os.path.exists(_install_dir):
            os.makedirs(_install_dir)
        _install_dir = os.path.abspath(_install_dir)

        # get output file path
        output = self.output()

        # set environment variables
        my_env = self.set_environment(self.env_script)

        # checking cmake path
        code, _cmake_executable, error = interruptable_popen(
            ["which", "cmake"], stdout=PIPE, stderr=PIPE, env=my_env
        )

        # actual payload:
        print("=========================================================")
        print("| Starting cmake step for CROWN")
        print("| Using cmake {}".format(_cmake_executable.replace("\n", "")))
        print("| Using CROWN {}".format(_crown_path))
        print("| Using build_directory {}".format(_build_dir))
        print("| Using install directory {}".format(_install_dir))
        print("=========================================================")

        # run CROWN build step
        _cmake_cmd = ["cmake", _crown_path]

        _cmake_args = [
            "-DANALYSIS={ANALYSIS}".format(ANALYSIS=_analysis),
            "-DSAMPLES={SAMPLES}".format(SAMPLES=_sampletype),
            "-DERAS={ERAS}".format(ERAS=_era),
            "-DCHANNELS={CHANNELS}".format(CHANNELS=_channels),
            "-DSHIFTS={SHIFTS}".format(SHIFTS=_shifts),
            "-DINSTALLDIR={INSTALLDIR}".format(INSTALLDIR=_install_dir),
            "-B{BUILDFOLDER}".format(BUILDFOLDER=_build_dir),
        ]
        print("Executable: {}".format(" ".join(_cmake_cmd + _cmake_args)))

        code, out, error = interruptable_popen(
            _cmake_cmd + _cmake_args, stdout=PIPE, stderr=PIPE, env=my_env
        )
        print(code, out, error)
        # if successful save Herwig-cache and run-file as tar.gz
        if code != 0:
            print("Error when running cmake {}".format(error))
            print("Output: {}".format(out))
            print("cmake returned non-zero exit status {}".format(code))
            raise Exception("cmake failed")
        else:
            print("Successful cmake build !")

        print("Executable: {}".format(" ".join(["make", "install"])))
        code, out, error = interruptable_popen(
            ["make", "install"],
            stdout=PIPE,
            stderr=PIPE,
            env=my_env,
            cwd=_build_dir,
        )
        if code != 0:
            print("Error when running make {}".format(error))
            print("Output: {}".format(out))
            print("make returned non-zero exit status {}".format(code))
            raise Exception("make failed")
        else:
            print("Successful cmake build !")

        # TODO Create Tarball from the install directory\
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
        print("Executable: {}".format(" ".join(command)))
        code, out, error = interruptable_popen(
            command,
            stdout=PIPE,
            stderr=PIPE,
            env=my_env,
            cwd=os.path.join(_install_dir),
        )
        if code != 0:
            print("Error when creating tarball {}".format(error))
            print("Output: {}".format(out))
            print("tar returned non-zero exit status {}".format(code))
            raise Exception("tar failed")
        else:
            print("Successful tarball creation !")
        # todo after done, move the created tarball into the tarballs folder

        print("=======================================================")
