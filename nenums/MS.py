#! /usr/bin/python3.6
# -*- coding: utf-8 -*-

"""
Class to work on a Measurement SET
        by A. Loh
"""

import os

from astropy.io import fits
from astropy.time import Time

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
        antTable(msname=self.msname, miniarrays=self.xst.miniarrays)

        emptyMS(msname=self.msname,
            start=Time(self.xst.xsttime[0], format='jd'),
            dt=self.xst.dt,
            bwidth=self.xst.bwidth,
            xstsbbands=self.xst.xstsubband )

        addInfos(msname=self.msname,
            xstheader=self.xst.head)

        addFreq(msname=self.msname, xstsbbands=self.xst.xstsubband )

        addTime(msname=self.msname, xsttime=self.xst.xsttime )

        addDescId(msname=self.msname, xstsbbands=self.xst.xstsubband )

        addData(msname=self.msname, builtma=self.xst.allmas, xstdata=self.xst.xstdata)

        zenithUVW(msname=self.msname)

        rephaseData(msname=self.msname, xsttime=self.xst.xsttime, ra_center=self.xst.ra, dec_center=self.xst.dec)

        addPointing(msname=self.msname, ra_center=self.xst.ra, dec_center=self.xst.dec)

        cleanDir(msname=self.msname)

        splitMS(msname=self.msname, remove=True)

        return


    # ================================================================= #
    # =========================== Internal ============================ #









