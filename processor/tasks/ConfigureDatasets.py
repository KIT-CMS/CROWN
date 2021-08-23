import law
import luigi
import os
import json
import yaml
from subprocess import PIPE
from law.util import interruptable_popen
from framework import Task


def get_energy(query):
    if "13TeV" in query:
        return 13
    elif "14TeV" in query:
        return 14


def get_era(query):
    if "UL16" in query:
        return "UL16"
    elif "UL17" in query:
        return "UL17"
    elif "UL18" in query:
        return "UL18"


def get_generator(query):
    pos = query.find("TeV") + 4
    generators = []
    if "amcatnlo" in query[pos:]:
        generators.append("amcatnlo")
    if "powheg" in query[pos:]:
        generators.append("powheg")
    if "madgraph" in query[pos:]:
        generators.append("madgraph")
    if "pythia" in query[pos:]:
        generators.append(
            query[pos:][query[pos:].find("pythia") : query[pos:].find("pythia") + 7]
        )
    generator = "-".join(generators)
    if generator == "":
        generator = "unspecified"
    return generator


def get_extension(query):
    startpos = query.find("ext")
    endpos = query.find("-v")
    return query[startpos:endpos]


def generate_nickname(query):
    subparts = query.split("/")
    nickname = subparts[1] + "-" + subparts[2]
    print(nickname)
    return nickname


def get_sample_type(datasetname):
    if "WJets" in datasetname:
        return "wjets"
    elif "Embedding" in datasetname:
        return "embedding"
    else:
        return "mc"


class ConfigureDatasets(Task):
    """
    Gather information on the selected datasets, for now just mentioning an explicit dataset
    """
    nick = luigi.Parameter()
    dataset = luigi.Parameter()
    env_script = os.path.join(
        os.path.dirname(__file__), "../../", "setup", "dasclient.sh"
    )

    def get_sample_information_from_das(self, nick, phys03):
        print("Getting sample information for \n  Nick: {}".format(nick))
        filedict = {}
        das_query = "dataset={}".format(self.dataset)
        if phys03:
            das_query += " instance=prod/phys03"
        else:
            das_query += " instance=prod/global"
        print("  DAS Query: {}".format(das_query))
        print("  Using Certificate {}".format(self.my_env["X509_USER_PROXY"]))
        cmd = [
            "/cvmfs/cms.cern.ch/common/dasgoclient --query '{}' --json".format(
                das_query
            )
        ]
        print("  Command: {}".format(" ".join(cmd)))
        code, out, error = interruptable_popen(
            cmd, stdout=PIPE, stderr=PIPE, shell=True, env=self.my_env
        )
        if code != 0:
            print("Error when trying to get DAS information {}".format(error))
            print("Output: {}".format(out))
            print("DAS query returned non-zero exit status {}".format(code))
            raise Exception("DAS query failed failed")
        result = json.loads(out)
        result = result[0]
        sample_dict = {}
        sample_dict["dbs"] = result["dataset"][0]["name"]
        sample_dict["nick"] = generate_nickname(sample_dict["dbs"])
        sample_dict["status"] = result["dataset"][0]["status"]
        sample_dict["datatier"] = result["dataset"][0]["data_tier_name"]
        sample_dict["prepid"] = result["dataset"][0]["prep_id"]
        sample_dict["datasetname"] = result["dataset"][0]["primary_dataset.name"]
        sample_dict["campaign"] = result["dataset"][0]["acquisition_era_name"]
        sample_dict["energy"] = get_energy(sample_dict["dbs"])
        sample_dict["era"] = get_era(sample_dict["dbs"])
        sample_dict["generators"] = get_era(sample_dict["dbs"])
        sample_dict["extension"] = get_extension(sample_dict["dbs"])
        sample_dict["version"] = result["dataset"][0]["processing_version"]
        sample_dict["sample_type"] = get_sample_type(sample_dict["datasetname"])

        return sample_dict

    def read_filelist_from_das(self, nick, phys03, xootd_prefix):
        filedict = {}
        das_query = "file dataset={}".format(self.dataset)
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
        filelist = json.loads(out)
        for file in filelist:
            filedict[file["file"][0]["name"]] = file["file"][0]["nevents"]
        print(
            "  Total files:  {} \n  Total events: {}".format(
                len(filedict.keys()), sum(filedict.values())
            )
        )
        return [
            "{prefix}/{path}".format(prefix=xootd_prefix, path=file)
            for file in filedict.keys()
        ], sum(filedict.values()), len(filedict.keys())

    def output(self):
        return self.local_target("sample_database/{}.yaml".format(self.nick))

    def run(self):
        # set environment variables
        self.my_env = self.set_environment(self.env_script)

        xootd_prefix = "root://cms-xrd-global.cern.ch/"
        xootd_prefix_gridka = "root://cmsxrootd-kit.gridka.de:1094/"
        """
         todo
        first check, if sample exists in datasets.yaml
            if it does, check if the corresponding filelist already exists, and if so, load this information as sample_data
            if not, query DAS, generate the filelist file and add all information to the datasets.yaml file

        Error handling for failing DAS query

        """
        output = self.output()
        output.parent.touch()
        # get general sample information from DAS
        sample_configuration = "sample_database/" + output.basename
        # if the filelist is already there, load it
        if os.path.exists(sample_configuration):
            print("Loading sample information from {}".format(sample_configuration))
            with open(sample_configuration, 'r') as stream:
                try:
                    sample_data = yaml.safe_load(stream)
                except yaml.YAMLError as exc:
                    print(exc)
                    raise Exception("Failed to load sample information")
        # otherwise, query DAS and generate the filelist
        else:
            sample_data = self.get_sample_information_from_das(self.dataset, phys03=False)
            sample_data["dbs"] = self.dataset
            sample_data["nick"] = self.dataset
            # read filelist information from DAS
            sample_data["filelist"], sample_data["nevents"], sample_data["nfiles"]= self.read_filelist_from_das(
                sample_data["nick"], False, xootd_prefix
            )
            # write the output

            file = open("sample_database/" + output.basename, "w")
            yaml.dump(sample_data, file)
            file.close()
        output.dump(sample_data)
