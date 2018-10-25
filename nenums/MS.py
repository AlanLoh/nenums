#! /usr/bin/python3.6
# -*- coding: utf-8 -*-

"""
Class to work on a Measurement SET
        by A. Loh
"""

import os

from astropy.io import fits

from .XST import XST
from .utils.ms import *

__author__ = ['Alan Loh']
__copyright__ = 'Copyright 2018, nenums'
__credits__ = ['Alan Loh']
__license__ = 'MIT'
__version__ = '0.0.1'
__maintainer__ = 'Alan Loh'
__email__ = 'alan.loh@obspm.fr'
__status__ = 'WIP'
__all__ = ['MS']


class MS(object):
    def __init__(self, xst, msname):
        self.xst    = xst
        self.msname = msname
    
    # ================================================================= #
    # ======================== Getter / Setter ======================== #
    @property
    def xst(self):
        """ XST observation
        """
        return self._xst
    @xst.setter
    def xst(self, x):
        self._xst = XST(xstfile=x)
        return


    # ================================================================= #
    # =========================== Methods ============================= #
    def createMS(self):
        """
        """
        antTable(msname=self.msname, miniarrays=self.miniarrays)

        emptyMS(msname=self.msname,
            start='2018-10-20 12:00:00',
            dt=1,
            bwidth=195.e3,
            xstsbbands=self.xstsubband )

        addInfos(msname=self.msname,
            xstheader=h)

        addFreq(msname=self.msname, xstsbbands=self.xstsubband )

        addTime(msname=self.msname, xsttime=self.xsttime )

        addDescId(msname=self.msname, xstsbbands=self.xstsubband )

        addData(msname=self.msname, builtma=self.allmas, xstdata=self.xstdata)

        zenithUVW(msname=self.msname)

        rephaseData(msname=self.msname, xsttime=self.xsttime, ra_center=45, dec_center=45)

        addPointing(msname=self.msname, ra_center=45, dec_center=45)

        cleanDir(msname=self.msname)

        splitMS(msname=self.msname, remove=True)

        return


    # ================================================================= #
    # =========================== Internal ============================ #









