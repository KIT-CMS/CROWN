# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

import os
import re
import json


class DatasetsHelper(object):
    def __init__(self, database_path):
        if not os.path.exists(database_path):
            logger.fatal(
                "Database file does not exist: {}".format(database_path))
            raise Exception
        from Kappa.Skimming.datasetsHelperTwopz import datasetsHelperTwopz
        self._datasets_helper = datasetsHelperTwopz(database_path)

    def get_nicks_with_query(self, query):
        return self._datasets_helper.get_nicks_with_query(query)


class DatasetsHelperLight(object):
    def __init__(self, database_path):
        self._database_path = database_path

    def _load_database(self):
        if not os.path.exists(self._database_path):
            logger.fatal("Database file does not exist: {}".format(
                self._database_path))
            raise Exception
        return json.load(open(self._database_path, "r"))

    def get_nicks_with_query(self, query):
        database = self._load_database()
        nicks = []
        for entry in database:
            passed = self._check_recursively(entry, query, database)
            if passed:
                nicks.append(entry)
        return nicks

    def _check_recursively(self, entry, query, database):
        for attribute in query:
            a = query[attribute]
            b = database[entry][attribute]
            if isinstance(b, str):
                result = re.match(a, b)
                if result == None:
                    return False
            elif isinstance(b, bool):
                if not a == b:
                    return False
            else:
                logger.critical("Can not process query with entry %s and %s.",
                                a, b)
                raise Exception
        return True
