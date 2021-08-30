import law
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
    dataset_database = luigi.Parameter()
    production_tag = luigi.Parameter()

    def requires(self):
        # load the list of samples to be processed
        with open(self.sample_list) as file:
            samples = [nick.replace("\n", "") for nick in file.readlines()]
        for nick in samples:
            # check if sample exists in datasets.yaml
            with open(self.dataset_database, "r") as stream:
                sample_db = yaml.safe_load(stream)
            if nick not in sample_db:
                console.log(
                    "Sample {} not found in {}".format(nick, self.dataset_database)
                )
                raise Exception("Sample not found in DB")
            sample_data = sample_db[nick]
            era = str(sample_data["era"])
            sampletype = sample_data["sample_type"]
            yield CROWNRun.req(self, nick=nick, era=era, sampletype=sampletype)

    def run(self):
        pass
