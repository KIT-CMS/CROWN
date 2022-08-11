import luigi
import yaml
from CROWNRun import CROWNRun
from framework import Task
from framework import console


class ProduceSamples(Task):
    """
    collective task to trigger ntuple production for a list of samples
    """

    sample_list = luigi.Parameter()
    analysis = luigi.Parameter()
    config = luigi.Parameter()
    dataset_database = luigi.Parameter()
    production_tag = luigi.Parameter()

    def requires(self):
        # load the list of samples to be processed
        data = {}
        data["sampletypes"] = set()
        data["eras"] = set()
        data["details"] = {}
        samples = []
        # check if sample list is a file or a comma separated list
        if self.sample_list.endswith(".txt"):
            with open(self.sample_list) as file:
                samples = [nick.replace("\n", "") for nick in file.readlines()]
        elif "," in self.sample_list:
            samples = self.sample_list.split(",")
        else:
            raise ValueError("sample-list is neither a file nor a comma separated list")
        for nick in samples:
            data["details"][nick] = {}
            # check if sample exists in datasets.yaml
            with open(self.dataset_database, "r") as stream:
                sample_db = yaml.safe_load(stream)
            if nick not in sample_db:
                console.log(
                    "Sample {} not found in {}".format(nick, self.dataset_database)
                )
                raise Exception("Sample not found in DB")
            sample_data = sample_db[nick]
            data["details"][nick]["era"] = str(sample_data["era"])
            data["details"][nick]["sampletype"] = sample_data["sample_type"]
            # all samplestypes and eras are added to a list, used to built the CROWN executable
            data["eras"].add(data["details"][nick]["era"])
            data["sampletypes"].add(data["details"][nick]["sampletype"])
        for nick in data["details"]:
            yield CROWNRun.req(
                self,
                nick=nick,
                era=data["details"][nick]["era"],
                sampletype=data["details"][nick]["sampletype"],
                eras=data["eras"],
                sampletypes=data["sampletypes"],
            )

    def run(self):
        pass
