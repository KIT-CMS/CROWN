# -*- coding: utf-8 -*-

import numpy as np

import logging
logger = logging.getLogger(__name__)
"""
"""


class Binning(object):
    pass


class VariableBinning(Binning):
    def __init__(self, bin_borders):
        # Convert bin borders to list of floats
        bin_borders = [float(f) for f in bin_borders]

        if not sorted(bin_borders) == bin_borders:
            logger.fatal("Binning %s is not sorted.", bin_borders)
            raise Exception
        if len(bin_borders) < 2:
            logger.fatal(
                "The binning has to consist of at least two border values.")
            raise Exception
        self._bin_borders = bin_borders

    @property
    def bin_borders(self):
        return np.array(self._bin_borders)

    @property
    def nbinsx(self):
        return len(self._bin_borders) - 1


class ConstantBinning(Binning):
    def __init__(self, nbinsx, xlow, xhigh):
        self._nbinsx = int(nbinsx)
        self._xlow = float(xlow)
        self._xhigh = float(xhigh)

    @property
    def bin_borders(self):
        return np.linspace(self._xlow, self._xhigh, self._nbinsx + 1)

    @property
    def nbinsx(self):
        return self._nbinsx

    @property
    def xlow(self):
        return self._xlow

    @property
    def xhigh(self):
        return self._xhigh
