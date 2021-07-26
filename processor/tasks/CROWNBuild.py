import law
import luigi
import os

from subprocess import PIPE
from law.util import interruptable_popen

from processor.framework import LocalTask, RemoteTask


class CROWNBuild(LocalTask):
    """
    Gather and compile CROWN with the given configuration
    """

    # configuration variables
    samples = luigi.Parameter()
    eras = luigi.Parameter()
    channels = luigi.Parameter()
    analysis = luigi.Parameter()
    build_cores = luigi.Parameter()
    shifts = luigi.Parameter()
    build_dir = luigi.Parameter()

    def convert_env_to_dict(self, env):
        my_env = {}
        for line in env.splitlines():
            if line.find(" ") < 0:
                try:
                    key, value = line.split("=", 1)
                    my_env[key] = value
                except ValueError:
                    pass
        return my_env

    def set_environment_variables(self):
        code, out, error = interruptable_popen(
            "source {}; env".format(
                os.path.join(
                    os.path.dirname(__file__),
                    "../../",
                    "setup",
                    "setup_crown_cmake.sh",
                )
            ),
            shell=True,
            stdout=PIPE,
            stderr=PIPE,
        )
        my_env = self.convert_env_to_dict(out)
        return my_env

    def output(self):
        return self.local_target("tarball.tar.gz")

    def create_build_directory(dir):
        dir = os.path.abspath(dir)
        if not os.path.exists(dir):
            os.makedirs(dir)

    def run(self):
        _samples = str(self.samples)
        _eras = str(self.eras)
        _channels = str(self.channels)
        _analysis = str(self.analysis)
        _build_cores = int(self.build_cores)
        _shifts = str(self.shifts)
        _build_dir = str(self.build_dir)

        # find crown
        _crown_path = os.path.abspath("CROWN")

        # create build directory
        if not os.path.exists(_build_dir):
            os.makedirs(_build_dir)
        _build_dir = os.path.abspath(_build_dir)

        # ensure that the output directory exists
        output = self.output()
        output.parent.touch()

        # set environment variables
        my_env = self.set_environment_variables()

        # checking cmake path
        code, _cmake_executable, error = interruptable_popen(
            ["which", "cmake"], stdout=PIPE, stderr=PIPE, env=my_env
        )

        # actual payload:
        print("=========================================================")
        print("| Starting cmake step for CROWN")
        print("| Using cmake {}".format(_cmake_executable))
        print("| Using CROWN {}".format(_crown_path))
        print("| Using build_directory {}".format(_build_dir))
        print("=========================================================")

        # run CROWN build step
        # _cd_builddir = ["cd" ,"{}".format(_build_dir), "&&"]
        # cmake command
        _cmake_cmd = ["cmake", _crown_path]

        _cmake_args = [
            "-DANALYSIS={ANALYSIS}".format(ANALYSIS=_analysis),
            "-DSAMPLES={SAMPLES}".format(SAMPLES=_samples),
            "-DERAS={ERAS}".format(ERAS=_eras),
            "-DCHANNELS={CHANNELS}".format(CHANNELS=_channels),
            "-DSHIFTS={SHIFTS}".format(SHIFTS=_shifts),
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

        print("Starting make with {} cores".format(_build_cores))
        print("Executable: {}".format(" ".join(["make", "install", "-j{}".format(_build_cores)])))
        code, out, error = interruptable_popen(
            ["make", "install", "-j{}".format(_build_cores)],
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

        print("=======================================================")
