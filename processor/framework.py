import os
import re
import luigi
import law
from law.util import interruptable_popen
from subprocess import PIPE
from rich.console import Console
from law.util import merge_dicts
import socket

law.contrib.load("wlcg")
law.contrib.load("htcondor")
console = Console(width=120)


class Task(law.Task):
    wlcg_path = luigi.Parameter()

    def store_parts(self):
        return (self.__class__.__name__,)

    def local_path(self, *path):
        parts = (os.getenv("ANALYSIS_DATA_PATH"),) + (self.__class__.__name__,) + path
        return os.path.join(*parts)

    def local_target(self, *path):
        return law.LocalFileTarget(self.local_path(*path))

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

    def set_environment(self, sourcescript):
        code, out, error = interruptable_popen(
            "source {}; env".format(sourcescript),
            shell=True,
            stdout=PIPE,
            stderr=PIPE,
        )
        if code != 0:
            console.log("Error when running source {}".format(error))
            console.log("Output: {}".format(out))
            console.log("source returned non-zero exit status {}".format(code))
            console.log("Script: {}".format(sourcescript))
            console.rule()
            raise Exception("source failed")
        my_env = self.convert_env_to_dict(out)
        return my_env

    def remote_path(self, *path):
        parts = (self.__class__.__name__,) + path
        return os.path.join(*parts)

    def remote_target(self, *path):
        target = law.wlcg.WLCGFileTarget(path=self.remote_path(*path))
        return target

    def remote_targets(self, paths):
        targets = []
        for path in paths:
            targets.append(law.wlcg.WLCGFileTarget(path=self.remote_path(path)))
        return targets


class HTCondorJobManager(law.htcondor.HTCondorJobManager):

    status_line_cre = re.compile(
        "^(\d+\.\d+)" + 4 * "\s+[^\s]+" + "\s+([UIRXSCHE<>])\s+.*$"
    )

    # def get_htcondor_version(cls):
    # return (9, 6, 0)
    #     return (8, 6, 5)


class HTCondorWorkflow(law.htcondor.HTCondorWorkflow):

    ENV_NAME = luigi.Parameter()
    htcondor_accounting_group = luigi.Parameter()
    htcondor_requirements = luigi.Parameter()
    htcondor_remote_job = luigi.Parameter()
    htcondor_user_proxy = luigi.Parameter()
    htcondor_walltime = luigi.Parameter()
    htcondor_request_cpus = luigi.Parameter()
    htcondor_request_gpus = luigi.Parameter(default="0")
    htcondor_request_memory = luigi.Parameter()
    htcondor_universe = luigi.Parameter()
    htcondor_docker_image = luigi.Parameter()
    htcondor_request_disk = luigi.Parameter()
    wlcg_path = luigi.Parameter()
    bootstrap_file = luigi.Parameter()
    replace_processor_tar = luigi.BoolParameter(default=False)

    # Use proxy file located in $X509_USER_PROXY or /tmp/x509up_u$(id) if empty
    htcondor_user_proxy = law.wlcg.get_voms_proxy_file()

    # output_collection_cls = law.SiblingFileCollection

    def htcondor_create_job_manager(self, **kwargs):
        kwargs = merge_dicts(self.htcondor_job_manager_defaults, kwargs)
        return HTCondorJobManager(**kwargs)

    def htcondor_output_directory(self):
        # Expand path to account for use of env variables (like $USER)
        return law.wlcg.WLCGDirectoryTarget(
            self.remote_path(),
            law.wlcg.WLCGFileSystem(
                None, base="{}".format(os.path.expandvars(self.wlcg_path))
            ),
        )

    def htcondor_create_job_file_factory(self):
        factory = super(HTCondorWorkflow, self).htcondor_create_job_file_factory()
        factory.is_tmp = False
        return factory

    def htcondor_bootstrap_file(self):
        hostfile = self.bootstrap_file
        return law.util.rel_path(__file__, hostfile)

    def htcondor_job_config(self, config, job_num, branches):
        analysis_name = os.getenv("ANA_NAME")
        # Check if env in config file is still the same as during the setup
        # TODO: is this assertion always intended?
        assert self.ENV_NAME == os.getenv(
            "ENV_NAME"
        ), "Environment of the config file ({}) is not the same as during the setup ({}).".format(
            self.ENV_NAME, os.getenv("ENV_NAME")
        )
        config.custom_content = []
        config.custom_content.append(
            ("accounting_group", self.htcondor_accounting_group)
        )
        config.render_variables["analysis_path"] = os.getenv("ANALYSIS_PATH")
        config.custom_content.append(("Requirements", self.htcondor_requirements))
        config.custom_content.append(("+RemoteJob", self.htcondor_remote_job))
        config.custom_content.append(("universe", self.htcondor_universe))
        config.custom_content.append(("docker_image", self.htcondor_docker_image))
        config.custom_content.append(("+RequestWalltime", self.htcondor_walltime))
        config.custom_content.append(("x509userproxy", self.htcondor_user_proxy))
        config.custom_content.append(("request_cpus", self.htcondor_request_cpus))
        # Only include "request_gpus" if any are requested, as nodes with GPU are otherwise excluded
        if float(self.htcondor_request_gpus) > 0:
            config.custom_content.append(("request_gpus", self.htcondor_request_gpus))
        config.custom_content.append(("RequestMemory", self.htcondor_request_memory))
        config.custom_content.append(("RequestDisk", self.htcondor_request_disk))

        prevdir = os.getcwd()
        os.system("cd $ANALYSIS_PATH")
        tarball = law.wlcg.WLCGFileTarget(
            path=self.remote_path(
                "job_tarball/{}_processor.tar.gz".format(analysis_name)
            )
        )
        if not os.path.exists("tarballs"):
            os.makedirs("tarballs")
        # TODO: how to determine if cfgs/tasks were changed and new tarball is necessary
        if (
            not os.path.isfile("tarballs/{}_processor.tar.gz".format(analysis_name))
            or self.replace_processor_tar
        ):
            command = [
                "tar",
                "--exclude",
                "*.pyc",
                "-czf",
                "tarballs/{}_processor.tar.gz".format(analysis_name),
                "processor",
                "lawluigi_configs/{}_luigi.cfg".format(analysis_name),
                "lawluigi_configs/{}_law.cfg".format(analysis_name),
                "law",
            ]
            code, out, error = interruptable_popen(
                command,
                stdout=PIPE,
                stderr=PIPE,
            )
            if code != 0:
                console.log("Error when taring job {}".format(error))
                console.log("Output: {}".format(out))
                console.log("tar returned non-zero exit status {}".format(code))
                console.rule()
                os.remove("tarballs/{}_processor.tar.gz".format(analysis_name))
                raise Exception("tar failed")
            else:
                console.rule("Successful tar!")
            console.rule("Uploading tarball to {}".format(tarball.path))
            tarball.parent.touch()
            tarball.copy_from_local(
                src="tarballs/{}_processor.tar.gz".format(analysis_name)
            )
            console.rule("Tarball uploaded!")
        # Check if env was found in cvmfs
        env_is_in_cvmfs = os.getenv("CVMFS_ENV_PRESENT")
        if env_is_in_cvmfs == "False":
            tarball_env = law.wlcg.WLCGFileTarget(
                path=self.remote_path("job_tarball/{}_env.tar.gz".format(self.ENV_NAME))
            )
            tarball_env.parent.touch()
            tarball_env.copy_from_local(
                src="tarballs/{}_env.tar.gz".format(self.ENV_NAME)
            )
        os.chdir(prevdir)
        # Make string with all environmental variables given to the job
        environment_string = ""
        environment_string += "ENV_NAME={} ".format(self.ENV_NAME)
        environment_string += "ANA_NAME={} ".format(os.getenv("ANA_NAME"))
        environment_string += "USER={} ".format(os.getenv("USER"))
        environment_string += "TARBALL_PATH={} ".format(
            os.path.expandvars(self.wlcg_path) + tarball.path
        )
        # Include path to env tarball if env not in cvmfs
        if env_is_in_cvmfs == "False":
            environment_string += "TARBALL_ENV_PATH={} ".format(
                os.path.expandvars(self.wlcg_path) + tarball_env.path
            )
        config.custom_content.append(("environment", '"' + environment_string + '"'))
        return config
