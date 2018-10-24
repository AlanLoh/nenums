#! /usr/bin/python3.6
# -*- coding: utf-8 -*-

"""
Class to work on a Measurement SET
        by A. Loh
"""

import os

from astropy.io import fits

from .XST import XST

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
    def __init__(self, xst):
        self.xst = xst
    
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

        hdu = fits.open('/data/loh/NenuFAR/Test_MS_labo/20180605_125000_XST.fits')
        s = hdu[7].data['xstsubband']
        t = hdu[7].data['jd']
        d = hdu[7].data['data']
        h = fits.getheader('/data/loh/NenuFAR/Test_MS_labo/20180605_125000_XST.fits', ext=0)
        m = np.squeeze(hdu[1].data['noMROn']) 
        allma = np.squeeze(hdu[1].data['noMR']) 

        antTable(msname='Test.ms', miniarrays=m)

        emptyMS(msname='Test.ms',
            start='2018-10-20 12:00:00',
            dt=1,
            bwidth=195.e3,
            xstsbbands=s)

        addInfos(msname='Test.ms',
            xstheader=h)

        addFreq(msname='Test.ms', xstsbbands=s)

        addTime(msname='Test.ms', xsttime=t)

        addDescId(msname='Test.ms', xstsbbands=s)

        addData(msname='Test.ms', builtma=allma, xstdata=d)

        zenithUVW(msname='Test.ms')

        #rephaseData(msname='Test.ms', xsttime=t, ra_center=45, dec_center=45)

        addPointing(msname='Test.ms', ra_center=45, dec_center=45)

        cleanDir(msname='Test.ms')

        splitMS(msname='Test.ms', remove=True)
        
        return


    # ================================================================= #
    # =========================== Internal ============================ #









