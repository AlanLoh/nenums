#! /usr/bin/python3.6
# -*- coding: utf-8 -*-

"""
Class work on a Measurement SET
        by A. Loh
"""

import os

from astropy.io import fits
from astropy.time import Time

try:
    from pyrap.tables import table, addImagingColumns
except:
    print("\n\t=== WARNING ===")
    # raise ImportError("\n\t=== Unable to import the pyrap module ===")


__author__ = ['Alan Loh']
__copyright__ = 'Copyright 2018, nenums'
__credits__ = ['Alan Loh']
__license__ = 'MIT'
__version__ = '0.0.1'
__maintainer__ = 'Alan Loh'
__email__ = 'alan.loh@obspm.fr'
__status__ = 'WIP'
__all__ = ['MS']


class MS():
    def __init__(self):
        return

    def emptyMS(self):
        """ Create an empty Measurement Set with the `makems` executable
            - Write a parset file 'makems.cfg': dictionnary of MS properties
            - Launch `makems`

        """
        if not os.path.isdir(self.msfile):
            #ra, dec = self._azel2radec(0, 90, Time(self.xsttime[0], format='jd')) # Zenith RA Dec at the beginning of the observation
            ra, dec = self._azel2radec(0, 90, self.xsttime[0]) # Zenith RA Dec at the beginning of the observation

            D = {}
            D['AntennaTableName']   = {'id': 0, 'val': os.path.join(self.savepath, 'ANTENNA')}
            D['Declination']        = {'id': 0, 'val': str( dec ) + 'rad'}
            D['RightAscension']     = {'id': 0, 'val': str( ra )  + 'rad'}
            D['MSName']             = {'id': 0, 'val': self.msfile}
            D['NBands']             = {'id': 0, 'val': str(self.subpertime )} # str(self.spec_windows.size)} # number of Spectral Windows
            D['WriteAutoCorr']      = {'id': 0, 'val': 'T'}
            D['StartFreq']          = {'id': 0, 'val': '[' + str( self.xstsubb[0, 0] * self.bandwidth + self.bandwidth/2. ) + ']'}
            D['StepFreq']           = {'id': 0, 'val': str(self.bandwidth)}
            #D['NFrequencies']       = {'id': 0, 'val': str(self.subpertime )} #str(self.spec_windows.size)} # number of distinct frequencies
            D['NFrequencies']       = {'id': 0, 'val': str(self.subbands.size )} #str(self.spec_windows.size)} # number of distinct frequencies
            D['StartTime']          = {'id': 0, 'val': Time(self.xsttime[0], format='jd').isot.replace('T', '/')}
            D['StepTime']           = {'id': 0, 'val': str(self.steptime)} #"."
            D['NTimes']             = {'id': 0, 'val': str(self.xstsubb.shape[0])}
            D['VDSPath']            = {'id': 0, 'val': self.savepath}
            D['WriteImagerColumns'] = {'id': 0, 'val': "T"}
            # ------ Create a parset from the dictionnary ------# 
            parset = open( os.path.join(self.savepath, 'makems.cfg'), 'w')
            for key in D.keys():
                parset.write('{} = {}\n'.format(key, D[key]['val']))
            parset.close()
            # ------ Use makems to produce an empty MS ------# 
            os.system('makems {}'.format(os.path.join(self.savepath, 'makems.cfg')))
            self._updateHistory(m='Creation of the Measurement Set from {}'.format(self.obsfile))
        else:
            print("\t=== Measurement Set {} already exists. ===".format(self.msfile))
        return

    def addInfos(self):
        """ 
        """
        return


    def addFreq(self):
        """
        """
        return




