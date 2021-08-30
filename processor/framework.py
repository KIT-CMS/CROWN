import os
import re
import luigi
import law
import law.contrib.htcondor
from law.util import interruptable_popen
from subprocess import PIPE
from rich.console import Console
from law.util import merge_dicts

law.contrib.load("wlcg")
console = Console(width=160)


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
        my_env = self.convert_env_to_dict(out)
        return my_env

    def remote_path(self, *path):
        parts = (self.__class__.__name__,) + path
        return os.path.join(*parts)

    def remote_target(self, *path):
        return law.wlcg.WLCGFileTarget(
            path=self.remote_path(*path),
            fs=law.wlcg.WLCGFileSystem(None, base="{}".format(self.wlcg_path)),
        )


class HTCondorJobManager(law.contrib.htcondor.HTCondorJobManager):

    status_line_cre = re.compile(
        "^(\d+\.\d+)" + 4 * "\s+[^\s]+" + "\s+([UIRXSCHE<>])\s+.*$"
    )

    def get_htcondor_version(cls):
        return (8, 6, 5)


class HTCondorWorkflow(law.contrib.htcondor.HTCondorWorkflow):

    htcondor_accounting_group = luigi.Parameter()
    htcondor_requirements = luigi.Parameter()
    htcondor_remote_job = luigi.Parameter()
    htcondor_user_proxy = luigi.Parameter()
    htcondor_walltime = luigi.Parameter()
    htcondor_request_cpus = luigi.Parameter()
    htcondor_request_memory = luigi.Parameter()
    htcondor_universe = luigi.Parameter()
    htcondor_docker_image = luigi.Parameter()
    htcondor_request_disk = luigi.Parameter()
    wlcg_path = luigi.Parameter()
    bootstrap_file_etp = luigi.Parameter()
    bootstrap_file_generic = luigi.Parameter()
    if "etp" in os.uname()[1]:
        etp = True
    else:
        etp = False

    output_collection_cls = law.SiblingFileCollection

    def htcondor_create_job_manager(self, **kwargs):
        kwargs = merge_dicts(self.htcondor_job_manager_defaults, kwargs)
        return HTCondorJobManager(**kwargs)

    def htcondor_output_directory(self):
        return law.wlcg.WLCGDirectoryTarget(
            self.remote_path(),
            law.wlcg.WLCGFileSystem(None, base="{}".format(self.wlcg_path)),
        )

    def htcondor_create_job_file_factory(self):
        factory = super(HTCondorWorkflow, self).htcondor_create_job_file_factory()
        factory.is_tmp = False
        return factory

    def htcondor_bootstrap_file(self):
        if self.etp:
            hostfile = self.bootstrap_file_etp
        else:
            hostfile = self.bootstrap_file_generic
        return law.util.rel_path(__file__, hostfile)

    def htcondor_job_config(self, config, job_num, branches):
        config.custom_content = []
        config.custom_content.append(
            ("accounting_group", self.htcondor_accounting_group)
        )
        # config.custom_content.append(("getenv", "true"))
        config.render_variables["analysis_path"] = os.getenv("ANALYSIS_PATH")
        config.custom_content.append(("Requirements", self.htcondor_requirements))
        config.custom_content.append(("+RemoteJob", self.htcondor_remote_job))
        config.custom_content.append(("universe", self.htcondor_universe))
        config.custom_content.append(("docker_image", self.htcondor_docker_image))
        config.custom_content.append(("+RequestWalltime", self.htcondor_walltime))
        config.custom_content.append(("x509userproxy", self.htcondor_user_proxy))
        config.custom_content.append(("request_cpus", self.htcondor_request_cpus))
        config.custom_content.append(("RequestMemory", self.htcondor_request_memory))
        config.custom_content.append(("RequestDisk", self.htcondor_request_disk))

        prevdir = os.getcwd()
        os.system("cd $ANALYSIS_PATH")
        if not self.etp:
            tarball = law.wlcg.WLCGFileTarget(
                path=self.remote_path("job_tarball/processor.tar.gz"),
                fs=law.wlcg.WLCGFileSystem(None, base="{}".format(self.wlcg_path)),
            )
        if not os.path.isfile("processor.tar.gz"):
            os.system(
                "tar --exclude *.pyc -czf processor.tar.gz processor luigi.cfg law.cfg law sample_database"
            )
            if not self.etp:
                tarball.parent.touch()
                tarball.copy_from_local(src="processor.tar.gz")
        os.chdir(prevdir)
        config.input_files.append(law.util.rel_path(__file__, "../processor.tar.gz"))
        if not self.etp:
            if not os.path.isfile("tarball_path.txt"):
                with open("tarball_path.txt", "w") as f:
                    f.write(self.wlcg_path + tarball.path)
            config.input_files.append(
                law.util.rel_path(__file__, "../tarball_path.txt")
            )
        return config
