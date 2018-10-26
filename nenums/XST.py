#! /usr/bin/python3.6
# -*- coding: utf-8 -*-

"""
Class to read and load XST data
        by A. Loh
"""

import os
import numpy as np
import glob

from astropy.io import fits
from astropy.time import Time

from .utils.astro import altaz2radec

__author__ = ['Alan Loh']
__copyright__ = 'Copyright 2018, nenums'
__credits__ = ['Alan Loh']
__license__ = 'MIT'
__version__ = '0.0.1'
__maintainer__ = 'Alan Loh'
__email__ = 'alan.loh@obspm.fr'
__status__ = 'WIP'
__all__ = ['XST']


class XST(object):
    def __init__(self, xstfile):
        self.xstfile = xstfile

    # ================================================================= #
    # ======================== Getter / Setter ======================== #
    @property
    def xstfile(self):
        """ XST observation file
        """
        return self._xstfile
    @xstfile.setter
    def xstfile(self, o):
        if (o is None) or (o == '') or (os.path.isdir(o)):
            if os.path.isdir(o):
                _opath = os.path.abspath(o)
                xstfiles = glob.glob( os.path.join(_opath, '*XST.fits') )
            else:
                xstfiles = glob.glob('*XST.fits')
            if len(xstfiles) == 0:
                raise IOError("\n\t=== No XST fits file in current directory, specify the file to read. ===")
            elif len(xstfiles) == 1:
                o = os.path.abspath(xstfiles[0])
            else:
                raise IOError("\n\t=== Multiple XST files are not handled yet ===")
        else:
            if not os.path.isfile(o):
                raise IOError("\t=== File {} not found. ===".format(os.path.abspath(o)))
            else:
                o = os.path.abspath(o) 
        self._xstfile = o

        if not self._isXST():
            raise ValueError("\t=== Files might not be XST observaiton ===")
        else:
            self._readXST()

        return

    # ================================================================= #
    # =========================== Internal ============================ #
    def _isXST(self):
        """ Check that self.xstfile is a proper XST observation
        """
        with fits.open(self.xstfile, mode='readonly', ignore_missing_end=True, memmap=True) as f:
            return f[0].header['OBJECT'] == 'crosscorrelation Statistics'

    def _readXST(self):
        """ Read XST fits files and fill the class attributes
        """
        self.obsname = os.path.basename( self.xstfile )[:-9]

        with fits.open(self.xstfile, mode='readonly', ignore_missing_end=True, memmap=True) as f:
            self.head = f[0].header
            setup_ins = f[1].data
            setup_ana = f[3].data
            setup_pbe = f[6].data
            data_hdu  = f[7].data

        self.miniarrays = np.array( [mrs[0:setup_ana['nbMRUsed'][i]] for i, mrs in enumerate(setup_ana['MRList'])] )[0]
        self.allmas     = setup_ins['noMR']
        self.xstsubband = data_hdu['xstsubband']
        self.xsttime    = data_hdu['jd']
        self.xstdata    = data_hdu['data']
        self.dt         = float(self.head['DT'])
        self.bwidth     = float(self.head['BANDWIDT'].split()[0])*1.e3 # in Hz
        self.ra, self.dec = altaz2radec(setup_pbe['AZ'][0], setup_pbe['EL'][0], Time(setup_pbe['timestamp'][0]))

        return



