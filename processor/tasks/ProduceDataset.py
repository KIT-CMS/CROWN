from KingMaker.processor.tasks.CROWNBuild import CROWNBuild
import law
import luigi
import os

from subprocess import PIPE
from law.util import interruptable_popen

from processor.framework import RemoteTask


class ProducerDataset(RemoteTask):
    """
    collective task to trigger ntuple production of a given dataset
    """

    dataset = luigi.Parameter()

    def output(self):
        return self.local_target("tarball.tar.gz")

    def requires(self):
        return CROWNBuild.req(self)

    def run(self):
        _dataset = str(self.dataset)

        # ensure that the output directory exists
        output = self.output()
        output.parent.touch()

        # set environment variables
        my_env = self.set_environment(self.env_script)

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

        print(
            "Executable: {}".format(
                " ".join(["make", "install", "-j{}".format(_build_cores)])
            )
        )
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
