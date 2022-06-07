import os
import re
import luigi
import law
from law.util import interruptable_popen
from subprocess import PIPE
from rich.console import Console
from law.util import merge_dicts
import socket
from datetime import datetime
from law.contrib.htcondor.job import HTCondorJobManager
from getpass import getuser

law.contrib.load("wlcg")
law.contrib.load("htcondor")
console = Console(width=120)

# Determine startup time to use as default production_tag
# LOCAL_TIMESTAMP is used by remote workflows to ensure consistent tags
if os.getenv('LOCAL_TIMESTAMP'):
    startup_time = os.getenv('LOCAL_TIMESTAMP')
else:
    startup_time = datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%f")

class Task(law.Task):

    wlcg_path = luigi.Parameter(
        description="Base-path to remote file location."
    )
    # Behaviour of production_tag:
    # If a tag is give it will be used for the respective task.
    # If no tag is given a timestamp abse on startup_time is used.
    #   This timestamp is the same for all tasks with no set production_tag.
    production_tag = luigi.Parameter(
        default="default/{}".format(startup_time), 
        description="Tag to differentiate workflow runs. Set to a timestamp as default."
    )
    output_collection_cls = law.NestedSiblingFileCollection

    # Path of local targets. Composed from the analysis path set during the setup.sh, 
    #   the production_tag, the name of the task and an additional path if provided.
    def local_path(self, *path):
        parts = (os.getenv("ANALYSIS_DATA_PATH"),) + (self.production_tag,) + (self.__class__.__name__,) + path
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

    # Function to apply a source-script and get the resulting environment.
    #   Anything apart from setting paths is likely not included in the resulting envs.
    def set_environment(self, sourcescript):
        console.log("Script: {}".format(sourcescript))
        if isinstance(sourcescript, str):
            sourcescript = [sourcescript]
        source_string = " ".join(sourcescript)
        code, out, error = interruptable_popen(
            "source {}; env".format(source_string),
            shell=True,
            stdout=PIPE,
            stderr=PIPE,
            # rich_console=console
        )
        if code != 0:
            console.log("source returned non-zero exit status {}".format(code))
            console.log("Error: {}".format(error))
            raise Exception("source failed")
        console.rule()
        my_env = self.convert_env_to_dict(out)
        return my_env

    # Run a bash command 
    #   Command can be composed of multiple parts (interpreted as seperated by a space).
    #   A sourcescript can be provided that is called by set_environment the resulting
    #       env is then used for the command
    #   The command is run as if it was called from run_location
    def run_command(self, command=[], sourcescript=[], run_location=None):
        if command:
            if isinstance(command, str):
                command = [command]
            if sourcescript:
                run_env = self.set_environment(sourcescript)
            else:
                run_env = None
            logstring = "Running {}".format(command)
            if run_location:
                logstring += " from {}".format(run_location)
            console.log(logstring)
            code, out, error = interruptable_popen(
                " ".join(command),
                shell=True,
                stdout=PIPE,
                stderr=PIPE,
                env=run_env,
                cwd=run_location,
                # rich_console=console
            )
            if code != 0:
                console.log("Error when running {}.".format(list(command)))
                console.log("Output: {}".format(out))
                console.log("Error: {}".format(error))
                console.log("Command returned non-zero exit status {}.".format(code))
                raise Exception("{} failed".format(list(command)))
            else:
                console.log("Command successful.")
                console.log("Output: {}".format(out))
                console.log("Error: {}".format(error))
        else:
            raise Exception("No command provided.")

    # Path of remote targets. Composed from the production_tag, 
    #   the name of the task and an additional path if provided.
    #   The wlcg_path will be prepended for WLCGFileTargets
    def remote_path(self, *path):
        parts = (self.production_tag,) + (self.__class__.__name__,) + path
        return os.path.join(*parts)

    def remote_target(self, *path):
        target = law.wlcg.WLCGFileTarget(path=self.remote_path(*path))
        return target

    def remote_targets(self, paths):
        targets = []
        for path in paths:
            targets.append(law.wlcg.WLCGFileTarget(path=self.remote_path(path)))
        return targets

class HTCondorWorkflow(Task, law.htcondor.HTCondorWorkflow):
    ENV_NAME = luigi.Parameter(
        description="Environment to be used in HTCondor job."
    )
    htcondor_accounting_group = luigi.Parameter(
        description="Accounting group to be set in Hthe TCondor job submission."
    )
    htcondor_requirements = luigi.Parameter(
        description="Job requirements to be set in the HTCondor job submission."
    )
    htcondor_remote_job = luigi.Parameter(
        description="Whether RemoteJob should be set in the HTCondor job submission."
    )
    htcondor_user_proxy = luigi.Parameter(
        description="VOMS-proxy in HTCondor job submission."
    )
    htcondor_walltime = luigi.Parameter(
        description="Runtime to be set in HTCondor job submission."
    )
    htcondor_request_cpus = luigi.Parameter(
        description="Number of CPU cores to be requested in HTCondor job submission."
    )
    htcondor_request_gpus = luigi.Parameter(
        default="0", 
        description="Number of GPUs to be requested in HTCondor job submission. Default is none."
    )
    htcondor_request_memory = luigi.Parameter(
        description="Amount of memory(MB) to be requested in HTCondor job submission."
    )
    htcondor_universe = luigi.Parameter(
        description="Universe to be set in HTCondor job submission."
    )
    htcondor_docker_image = luigi.Parameter(
        description="Docker image to be used in HTCondor job submission."
    )
    htcondor_request_disk = luigi.Parameter(
        description="Amount of scratch-space(kB) to be requested in HTCondor job submission."
    )
    bootstrap_file = luigi.Parameter(
        description="Bootstrap script to be used in HTCondor job to set up law."
    )
    additional_files = luigi.ListParameter(
        default=[], 
        description="Additional files to be included in the job tarball. Will be unpacked in the run directory"
    ) 

    # Use proxy file located in $X509_USER_PROXY or /tmp/x509up_u$(id) if empty
    htcondor_user_proxy = law.wlcg.get_voms_proxy_file()

    def htcondor_create_job_manager(self, **kwargs):
        kwargs = merge_dicts(self.htcondor_job_manager_defaults, kwargs)
        return HTCondorJobManager(**kwargs)

    def htcondor_output_directory(self):
        # Expand path to account for use of env variables (like $USER)
        return law.wlcg.WLCGDirectoryTarget(
            self.remote_path("htcondor_files"),
            law.wlcg.WLCGFileSystem(
                None, base="{}".format(os.path.expandvars(self.wlcg_path))
            ),
        )

    def htcondor_create_job_file_factory(self):
        factory = super(HTCondorWorkflow, self).htcondor_create_job_file_factory()
        factory.is_tmp = False
        # Print location of job dir
        console.log("HTCondor job directory is: {}".format(factory.dir))
        return factory

    def htcondor_bootstrap_file(self):
        hostfile = self.bootstrap_file
        return law.util.rel_path(__file__, hostfile)
        
    def htcondor_job_config(self, config, job_num, branches):
        analysis_name = os.getenv("ANA_NAME")
        task_name = self.__class__.__name__
        analysis_path = os.getenv("ANALYSIS_PATH")
        # Write job config file 
        config.custom_content = []
        config.custom_content.append(
            ("accounting_group", self.htcondor_accounting_group)
        )
        # config.custom_content.append(("Log", "log.txt")) #
        # config.custom_content.append(("stream_output", "True")) #
        # config.custom_content.append(("Output", "out_{}to{}.txt".format(branches[0], branches[-1]))) #Remove before commit
        # config.custom_content.append(("stream_error", "True")) #
        # config.custom_content.append(("Output", "err_{}to{}.txt".format(branches[0], branches[-1]))) #
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

        # Ensure tarball dir exists
        if not os.path.exists("tarballs/{}".format(self.production_tag)):
            os.makedirs("tarballs/{}".format(self.production_tag))
        # Repack tarball if it is not available remotely
        tarball = law.wlcg.WLCGFileTarget(
            "{tag}/{task}/job_tarball/processor.tar.gz".format(
                tag=self.production_tag,
                task=self.__class__.__name__
            )
        )
        if not tarball.exists():
            # Make new tarball 
            prevdir = os.getcwd()
            os.system("cd $ANALYSIS_PATH")
            tarball_local = law.LocalFileTarget(
                "tarballs/{}/{}/processor.tar.gz".format(
                    self.production_tag, task_name
                )
            )
            tarball_local.parent.touch()
            # Create tarball containing:
            #   The processor directory, thhe relevant config files, law
            #   and any other files specified in the additional_files parameter
            command = [
                "tar",
                "--exclude",
                "*.pyc",
                "--exclude",
                "law/.git",
                "-czf",
                "tarballs/{}/{}/processor.tar.gz".format(
                    self.production_tag, task_name
                ),
                "processor",
                "lawluigi_configs/{}_luigi.cfg".format(analysis_name),
                "lawluigi_configs/{}_law.cfg".format(analysis_name),
                "law",
            ] + list(self.additional_files)
            code, out, error = interruptable_popen(
                command,
                stdout=PIPE,
                stderr=PIPE,
                # rich_console=console
            )
            if code != 0:
                console.log("Error when taring job {}".format(error))
                console.log("Output: {}".format(out))
                console.log("tar returned non-zero exit status {}".format(code))
                console.rule()
                os.remove(
                    "tarballs/{}/{}/processor.tar.gz".format(
                        self.production_tag, task_name
                    )
                )
                raise Exception("tar failed")
            else:
                console.rule("Successful tar!")
            # Copy new tarball to remote
            tarball.parent.touch()
            tarball.copy_from_local(
                src="tarballs/{}/{}/processor.tar.gz".format(
                    self.production_tag, task_name
                )
            )
            console.rule("Tarball uploaded!")
            os.chdir(prevdir)
        # Check if env of this task was found in cvmfs
        env_list = os.getenv('ENV_NAMES_LIST').split(';')
        env_list = list(dict.fromkeys(env_list[:-1]))
        env_dict = dict(env.split(',') for env in env_list)
        if env_dict[self.ENV_NAME] == "False":
            # IMPORTANT: environments have to be named differently with each change
            #            as chaching prevents a clean overwrite of existing files
            tarball_env = law.wlcg.WLCGFileTarget(
                path="env_tarballs/{env}.tar.gz".format(
                    env=self.ENV_NAME
                )
            )
            if not tarball_env.exists():
                tarball_env.parent.touch()
                tarball_env.copy_from_local(
                    src="tarballs/conda_envs/{}.tar.gz".format(self.ENV_NAME)
                )
        config.render_variables["USER"] = getuser()
        config.render_variables["ANA_NAME"] = os.getenv("ANA_NAME")
        config.render_variables["ENV_NAME"] = self.ENV_NAME
        config.render_variables["TAG"] = self.production_tag
        config.render_variables["USE_CVMFS"] = env_dict[self.ENV_NAME]
        config.render_variables["TARBALL_PATH"] = os.path.expandvars(self.wlcg_path) + tarball.path
        # Include path to env tarball if env not in cvmfs
        if env_dict[self.ENV_NAME] == "False":
            config.render_variables["TARBALL_ENV_PATH"] = os.path.expandvars(self.wlcg_path) + tarball_env.path
        config.render_variables["LOCAL_TIMESTAMP"] = startup_time

        return config
