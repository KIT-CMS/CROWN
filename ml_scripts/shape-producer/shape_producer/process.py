# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)
"""
"""


class Process(object):
    def __init__(self, name, estimation):
        self._name = name
        self._estimation_method = estimation

    @property
    def name(self):
        return self._name

    @property
    def estimation_method(self):
        return self._estimation_method


class Processes(object):
    def __init__(self):
        self._processes = []

    @property
    def processes(self):
        return self._processes

    def select(self, *args):
        return [p for p in self._processes if p.name in args]

    def add(self, process):
        self._processes.append(process)
