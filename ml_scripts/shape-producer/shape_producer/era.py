# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

from .cutstring import Constant
from .datasets_helper import DatasetsHelperLight as DatasetsHelper

import os
"""
"""


def log_query(name, query, files):
    logger.debug("Query for %s: %s", name, query)
    if len(files) < 1:
        logger.critical("Query for %s did not return any files.", name)
        raise Exception
    for i_file, file_ in enumerate(files):
        logger.debug("File %d: %s", i_file + 1, file_)


class Era(object):
    def __init__(self, name, luminosity, database_path):
        self._name = name
        self._luminosity = luminosity
        self._database_path = database_path
        self._datasets_helper = DatasetsHelper(self._database_path)

    @property
    def datasets_helper(self):
        return self._datasets_helper

    @property
    def name(self):
        return self._name

    @property
    def lumi_weight(self):
        return Constant(str(self._luminosity), "lumi")


class Run2016(Era):
    def __init__(self, database_path):
        super(Run2016, self).__init__("Run2016", 35.87 * 1000.0, database_path)

    def data_files(self, channel):
        query = {
            "data": True,
            "campaign": "Run2016(B|C|D|E|F|G|H)",
            "scenario": "17Jul2018.*"
        }
        if channel.name == "mt" or channel.name == "mm":
            query["process"] = "SingleMuon"
        elif channel.name == "et":
            query["process"] = "SingleElectron"
        elif channel.name == "tt":
            query["process"] = "Tau"
        elif channel.name == "em":
            query["process"] = "MuonEG"
        elif channel.name == "ee":
            query["process"] = "DoubleEG"
        else:
            logger.critical("Channel %s is not implemented.", channel.name)
        files = self.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return files

class Run2017(Era):
    def __init__(self, database_path):
        super(Run2017, self).__init__("Run2017",
                                                 41.529 * 1000.0, database_path)

    def data_files(self, channel):
        query = {
            "data": True,
            "campaign": "Run2017(B|C|D|E|F)",
            "scenario": "31Mar2018v1"
        }
        if channel.name == "mt" or channel.name == "mm":
            query["process"] = "SingleMuon"
        elif channel.name == "et":
            query["process"] = "SingleElectron"
        elif channel.name == "tt":
            query["process"] = "Tau"
        elif channel.name == "em":
            query["process"] = "MuonEG"
        else:
            logger.critical("Channel %s is not implemented.", channel.name)
        files = self.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return files


class Run2018(Era):
    def __init__(self, database_path):
        super(Run2018, self).__init__("Run2018", 59.7 * 1000.0, database_path)

    def data_files(self, channel):
        query = {
            "data": True,
            "campaign": "Run2018(A|B|C|D)"
        }
        if channel.name == "mt":
            query["process"] = "SingleMuon"
        elif channel.name == "et":
            query["process"] = "EGamma"
        elif channel.name == "tt":
            query["process"] = "Tau"
        elif channel.name == "em":
            query["process"] = "MuonEG"
        elif channel.name == "mm":
            query["process"] = "SingleMuon"
        else:
            logger.critical("Channel %s is not implemented.", channel.name)
        files = self.datasets_helper.get_nicks_with_query(query)
        log_query(self.name, query, files)
        return files
