#! /usr/bin/python3.6
# -*- coding: utf-8 -*-

"""
Class to read and load XST data
        by A. Loh
"""

import os

from astropy.io import fits
from astropy.time import Time

from . import MS

__author__ = ['Alan Loh']
__copyright__ = 'Copyright 2018, nenums'
__credits__ = ['Alan Loh']
__license__ = 'MIT'
__version__ = '0.0.1'
__maintainer__ = 'Alan Loh'
__email__ = 'alan.loh@obspm.fr'
__status__ = 'WIP'
__all__ = ['XST']


class XST()
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
                bstfiles = glob.glob( os.path.join(_opath, '*XST.fits') )
            else:
                bstfiles = glob.glob('*XST.fits')
            if len(bstfiles) == 0:
                raise IOError("\n\t=== No XST fits file in current directory, specify the file to read. ===")
            elif len(bstfiles) == 1:
                o = os.path.abspath(bstfiles[0])
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
    # =========================== Converter =========================== #
    def convertMS(self):
        """
        """
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
            head      = f[0].header
            setup_ins = f[1].data
            setup_obs = f[2].data
            setup_ana = f[3].data
            setup_bea = f[4].data
            setup_pan = f[5].data
            setup_pbe = f[6].data

        self.obstart  = Time( head['DATE-OBS'] + 'T' + head['TIME-OBS'] )
        self.obstop   = Time( head['DATE-END'] + 'T' + head['TIME-END'] )
        self.time     = [self.obstart.copy(), self.obstop.copy()]
        self.exposure = self.obstop - self.obstart

        self.miniarrays = np.squeeze( setup_ins['noMROn'] )
        self._marot     = np.squeeze( setup_ins['rotation'] )
        self._mapos     = np.squeeze( setup_ins['noPosition'] )
        self._mapos     = self._mapos.reshape( int(self._mapos.size/3), 3 )
        self._pols      = np.squeeze( setup_ins['spol'] )

        self.abeams     = setup_ana['NoAnaBeam']
        self._antlist   = np.array( [ np.array(eval(i)) - 1 for i in setup_ana['Antlist'] ])

        self.dbeams     = setup_bea['noBeam']
        self._digi2ana  = setup_bea['NoAnaBeam']
        self._bletlist  = setup_bea['BeamletList']
        self._freqs     = setup_bea['freqList']

        self._pointana  = setup_pan['noAnaBeam']
        self._azlistana = setup_pan['AZ']
        self._ellistana = setup_pan['EL']
        self._pointanat = Time(np.array([setup_pan['timestamp'][self._pointana==i] for i in self.abeams]))
        
        self._pointdig  = setup_pbe['noBeam']
        self._azlistdig = setup_pbe['AZ']
        self._ellistdig = setup_pbe['EL']    
        self._pointdigt = Time(setup_pbe['timestamp'])
        return



