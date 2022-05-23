# -*- coding: utf-8 -*-

from cutstring import *
import logging
logger = logging.getLogger(__name__)
"""
"""


# Class to store a variable and its binning
class Variable(object):
    def __init__(self, name, binning, expression=None):
        self._name = name
        self._binning = binning
        if expression == None:
            self._expression = name
        else:
            self._expression = expression

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return super(type(self), self).__repr__().split(' object at')[0][1:] + '(%s, %s) [addr.: %s]' % (self._name, self._expression, hex(id(self)))
    @property
    def name(self):
        return self._name

    @property
    def binning(self):
        return self._binning

    @property
    def expression(self):
        return self._expression
