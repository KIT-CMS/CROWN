import luigi
import os
import json
import yaml
from subprocess import PIPE
from law.util import interruptable_popen
from framework import Task
from framework import console


def ensure_dir(file_path):
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


class ConfigureDatasets(Task):
    """
    Gather information on the selected datasets, for now just mentioning an explicit dataset
    """

    nick = luigi.Parameter()
    dataset_database = luigi.Parameter()
    production_tag = luigi.Parameter()
    env_script = os.path.join(
        os.path.dirname(__file__), "../../", "setup", "dasclient.sh"
    )

    def read_filelist_from_das(self, dbs, phys03, xootd_prefix):
        filedict = {}
        das_query = "file dataset={}".format(dbs)
        if phys03:
            das_query += " instance=prod/phys03"
        else:
            das_query += " instance=prod/global"
        cmd = [
            "/cvmfs/cms.cern.ch/common/dasgoclient --query '{}' --json".format(
                das_query
            )
        ]
        code, out, error = interruptable_popen(
            cmd, stdout=PIPE, stderr=PIPE, shell=True, env=self.my_env
        )
        # jsonS = code.communicate()[0]
        if code != 0:
            raise Exception(
                "DAS query failed: {} \n Error: {}".format(das_query, error)
            )
        filelist = json.loads(out)
        for file in filelist:
            filedict[file["file"][0]["name"]] = file["file"][0]["nevents"]
        return (
            [
                "{prefix}/{path}".format(prefix=xootd_prefix, path=file)
                for file in filedict.keys()
            ],
            sum(filedict.values()),
            len(filedict.keys()),
        )

    def output(self):
        target = self.remote_target("sample_database/{}.yaml".format(self.nick))
        return target

    def run(self):
        output = self.output()
        output.parent.touch()
        if not output.exists():
            # set environment variables
            self.my_env = self.set_environment(self.env_script)

            # in germany, the european is the fastest
            xootd_prefix_global = "root://cms-xrd-global.cern.ch/"
            xootd_prefix_gridka = "root://cmsxrootd-kit.gridka.de:1094/"
            xootd_prefix_europe = "root://xrootd-cms.infn.it/"

            output.parent.touch()

            with open(self.dataset_database, "r") as stream:
                sample_db = yaml.safe_load(stream)
            sample_data = sample_db[self.nick]

            sample_configfile = "sample_database/{era}/{type}/{nick}.yaml".format(
                era=sample_data["era"], type=sample_data["sample_type"], nick=self.nick
            )
            console.rule("Nick: {}".format(self.nick))
            # if the filelist is already there, load it
            if os.path.exists(sample_configfile):
                print("Loading sample config file")
                with open(sample_configfile, "r") as stream:
                    try:
                        sample_data = yaml.safe_load(stream)
                    except yaml.YAMLError as exc:
                        print(exc)
                        raise Exception("Failed to load sample information")
            # otherwise, query DAS and generate the filelist
            else:
                print("Loading from DAS")
                sample_data = sample_db[self.nick]
                sample_data["nick"] = self.nick
                # read filelist information from DAS
                (
                    sample_data["filelist"],
                    sample_data["nevents"],
                    sample_data["nfiles"],
                ) = self.read_filelist_from_das(
                    sample_data["dbs"], False, xootd_prefix_europe
                )
                # write the output
                ensure_dir(sample_configfile)
                file = open(sample_configfile, "w")
                yaml.dump(sample_data, file)
                file.close()

            console.log("Total Files: {}".format(sample_data["nfiles"]))
            console.log("Total Events: {}".format(sample_data["nevents"]))
            console.rule()
            output.dump(sample_data)
