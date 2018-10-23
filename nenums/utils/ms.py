#! /usr/bin/python3.6
# -*- coding: utf-8 -*-

"""
Measurement Set creation functions
        by A. Loh
"""

import os
import numpy as np
import warnings
import inspect

from astropy.time import Time
try:
    from pyrap.tables import table, addImagingColumns
except:
    print("\n\t=== WARNING ===")

from .astro import altaz2radec

__author__ = 'Alan Loh'
__copyright__ = 'Copyright 2018, nenums'
__credits__ = ['Alan Loh']
__version__ = '0.0.1'
__maintainer__ = 'Alan Loh'
__email__ = 'alan.loh@obspm.fr'
__status__ = 'Production'
__all__ = ['emptyMS',
            'addInfos',
            'addFreq',
            'addData',
            'zenithUVW',
            'rephaseData',
            'updateHist']


# =================================================================================== #
# ------------------------------------- emptyMS ------------------------------------- #
# =================================================================================== #
def emptyMS(msname, outpath, start, dt, bwidth, sbpertime, subbands):
    """ Create an empty *Measurement Set* designed for NenuFAR XST observations.

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **outpath** : str
            Path to write the Measurmeent Set
        * **start** : (str, Time)
            Beginning of the observation. XST cross-correlations are computed at the local zenith.
            It helps initializing the tracked RA / Dec.
        * **dt** : float
            Time step of the XST observation, in seconds
        * **bwidth** : float
            Spectral bandwidth in Hz
        * **sbpertime** : int
            Number of sub-bands per unit of time
    """
    if os.path.isdir(msname):
        warnings.warn('Measurement Set {} already exist'.format(msname))
        return
    if isinstance(start, str):
        try:
            start = Time(start)
        except ValueError:
            print("start syntax is maybe incorrect --> 'YYYY-MM-DD hh:mm:ss'")
    elif isinstance(start, Time):
        pass
    else:
        raise ValueError('Unrecognized time format for start keyword')

    ra, dec = altaz2radec(0., 90., start)

    parset_dict = {}
    parset_dict['AntennaTableName']   = {'id': 0, 'val': os.path.join(outpath, 'ANTENNA')}
    parset_dict['Declination']        = {'id': 0, 'val': str(dec) + 'rad'}
    parset_dict['RightAscension']     = {'id': 0, 'val': str(ra)  + 'rad'}
    parset_dict['MSName']             = {'id': 0, 'val': msname}
    parset_dict['NBands']             = {'id': 0, 'val': str(sbpertime)}
    parset_dict['WriteAutoCorr']      = {'id': 0, 'val': 'T'}
    parset_dict['StartFreq']          = {'id': 0, 'val': '[' + str(subbands[0, 0]*bwidth + bwidth/2.) + ']'}
    parset_dict['StepFreq']           = {'id': 0, 'val': str(bwidth)}
    parset_dict['NFrequencies']       = {'id': 0, 'val': str(subbands.size )}
    parset_dict['StartTime']          = {'id': 0, 'val': start.isot.replace('T', '/')}
    parset_dict['StepTime']           = {'id': 0, 'val': str(dt)}
    parset_dict['NTimes']             = {'id': 0, 'val': str(subbands.shape[0])}
    parset_dict['VDSPath']            = {'id': 0, 'val': outpath}
    parset_dict['WriteImagerColumns'] = {'id': 0, 'val': "T"}
    
    # ------ Create a parset from the dictionnary ------# 
    parset = open( os.path.join(outpath, 'makems.cfg'), 'w')
    for key in parset_dict.keys():
        parset.write('{} = {}\n'.format(key, parset_dict[key]['val']))
    parset.close()
    
    # ------ Use makems to produce an empty MS ------# 
    # os.system('makems {}'.format(os.path.join(self.savepath, 'makems.cfg')))
        
    # updateHist()

    return


# =================================================================================== #
# ------------------------------------ addInfos ------------------------------------- #
# =================================================================================== #
def addInfos(msname, xstheader):
    """ 
    """
    mstable   = table( os.path.join(msname, 'OBSERVATION'), ack=False, readonly=False)
    observer  = mstable.getcol('OBSERVER')
    project   = mstable.getcol('PROJECT')
    schedule  = mstable.getcol('SCHEDULE_TYPE')
    telescope = mstable.getcol('TELESCOPE_NAME')

    observer[0]  = xstheader['CONTACT']
    project[0]   = xstheader['OBJECT']
    schedule[0]  = xstheader['INSTRUME']
    telescope[0] = xstheader['INSTRUME']

    mstable.putcol('OBSERVER', observer)
    mstable.putcol('PROJECT', project)
    mstable.putcol('SCHEDULE_TYPE', schedule)
    mstable.putcol('TELESCOPE_NAME', telescope)

    mstable.flush()
    mstable.close()

    del observer, project, schedule, telescope, mstable

    updateHist(message='Observation info added, from XST fits (version {})\
        created by software version {}.'.format(xstheader['VERSION'], xstheader['SOFTWARE']))
    
    return


# =================================================================================== #
# ------------------------------------- addFreq ------------------------------------- #
# =================================================================================== #
def addFreq():
    """
    """
    return


# =================================================================================== #
# ------------------------------------- addTime ------------------------------------- #
# =================================================================================== #
def addTime(msname, xsttime, dt, sbpertime):
    """ Fill the TIME column of the MS based on XST times (in Julian Days)

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **xsttime** : np.ndarray
            Time array from a XST file (in JD)
        * **dt** : float
            Time step of the XST observation, in seconds
        * **sbpertime** : int
            Number of sub-bands per unit of time
    """
    _isMS(msname)

    if not isinstance(sbpertime, (int, np.integer)):
        raise TypeError("\t=== 'sbpertime' must be an integer ===")

    antable = table( os.path.join(msname, 'ANTENNA'), ack=False, readonly=True)
    baselines = int(antable.nrows() * (antable.nrows() - 1)/2. + antable.nrows())
    antable.close()

    mstable    = table(msname, ack=False, readonly=False)
    time_table = Time(xsttime, format='jd').mjd * 3600. * 24. + (dt/2.) # MJD in sec (MS standards)
    time_table = np.repeat(time_table, baselines * sbpertime)
    mstable.putcol('TIME', time_table)  
    mstable.flush()
    mstable.close()

    updateHist(message='TIME table filled.')
    return



# =================================================================================== #
# ------------------------------------ addDescId ------------------------------------ #
# =================================================================================== #
def addDescId(msname, xstsbbands):
    """ Fill the Data Description ID column of the MS based on XST subbands info

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **xstsbbands** : np.ndarray
            Sub-bands of a XST file (ext: 7, key: 'xstsubband')
    """
    _isMS(msname)

    subbands = np.unique(xstsbbands)

    antable = table( os.path.join(msname, 'ANTENNA'), ack=False, readonly=True)
    baselines = int(antable.nrows() * (antable.nrows() - 1)/2. + antable.nrows())
    antable.close()

    mstable = table(msname, ack=False, readonly=False)
    desc_table   = mstable.getcol('DATA_DESC_ID')
    spw_idx_time = np.searchsorted(subbands, xstsbbands.ravel()) # associate the SPW index to observed Sub-Band
    # Make sure that unused SPW are not indexed
    mask         = spw_idx_time < subbands.size
    mask[mask]   = (subbands[spw_idx_time[mask]] == xstsbbands.ravel()[mask])
    spw_idx_time = spw_idx_time[mask] # SPW index per unit of time
    spw_idx_time = np.repeat(spw_idx_time, baselines) # repeat for the number of visibilities
    mstable.putcol('DATA_DESC_ID', spw_idx_time)  
    mstable.flush()
    mstable.close()

    updateHist(message='TIME table filled.')
    return


# =================================================================================== #
# ------------------------------------- addData ------------------------------------- #
# =================================================================================== #
def addData(msname, builtma, xstdata):
    """ Fill the DATA column of the MS based on XST cross-correlations
        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **builtma** : np.ndarray
            Mini-Arrays on the field, from the XST file (ext: 1, key: 'noMR')
        * **xstdata** : np.ndarray
            Cross-correlations from the XST file (ext: 7, key: 'data')
    """
    antable   = table( os.path.join(msname, 'ANTENNA'), ack=False, readonly=True)
    baselines = int(antable.nrows() * (antable.nrows() - 1)/2. + antable.nrows())
    antnames  = antable.getcol('NAME')
    usedmas   = np.array([ int(mr.replace('MR', '')) for mr in antnames ])
    antable.close()

    mstable = table(msname, ack=False, readonly=False)
    data_ms = mstable.getcol('DATA')
    ant1_table = mstable.getcol('ANTENNA1')
    ant2_table = mstable.getcol('ANTENNA2')

    # bar = ChargingBar('Filling MS data from the XST cross-correlations', max=baselines)
    for ant1 in builtma:
        for ant2 in builtma[ant1:]:
            if (ant1 not in usedmas) or (ant2 not in usedmas):
                continue
            else:
                # Get the time x frequency x correl for the current baseline
                baseline = np.sort( [ant1, ant2] )[::-1] # descending order
                # These formulas were defined to fit the table of detail_FITS_00_05.pdf (kinda triangular series)
                MA1X_MA2X = baseline[0]*2*(baseline[0]*2+1)/2       + 2*baseline[1]
                MA1X_MA2Y = baseline[0]*2*(baseline[0]*2+1)/2+1     + 2*baseline[1] # it will be the same as next one for auto-corr
                MA1Y_MA2X = (baseline[0]*2+1)*(baseline[0]*2+2)/2   + 2*baseline[1]
                MA1Y_MA2Y = (baseline[0]*2+1)*(baseline[0]*2+2)/2+1 + 2*baseline[1]
                index     = np.array([MA1X_MA2X, MA1X_MA2Y, MA1Y_MA2X, MA1Y_MA2Y])
                data      = xstdata[:, :, index] # (time, frequency, 4)
                if ant1 == ant2:
                    # Index 1 (XY) is wrong for auto-correlations with the above algorithm,
                    # doesn't matter, we replace the data by the conjugate of YX which exists anyway.
                    data[:, :, 1] = np.conjugate(data[:, :, 2]) 
                ant1_index   = np.where(usedmas == ant1)[0][0]
                ant2_index   = np.where(usedmas == ant2)[0][0]
                clever_index = (ant1_table == ant1_index) & (ant2_table == ant2_index)
                data_ms[clever_index, :, :] = np.reshape(data, (data.shape[0]*data.shape[1], 1, data.shape[2]))
                # bar.next()
    # bar.finish()

    mstable.putcol('DATA', data_ms) 
    mstable.flush()
    mstable.close()

    _checkSum()

    updateHist(message='DATA table filled with XST cross-correlations.')
    return


# =================================================================================== #
# ----------------------------------- zenithUVW ------------------------------------- #
# =================================================================================== #
def zenithUVW(msname, sbpertime):
    """ Re-compute the UVW table of the Measurement Set.
        `emptyMS()` automatically creates a UVW table as if initial RA Dec (zenith at start) was tracked.
        However, XST are cross-correlations fixed at the local zenith.
        Therefore, the correct UVW table must be copied from the first time step.

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **sbpertime** : int
            Number of sub-bands per unit of time
    """
    _isMS(msname)

    if not isinstance(sbpertime, (int, np.integer)):
        raise TypeError("\t=== 'sbpertime' must be an integer ===")
    
    mstable = table(msname, ack=False, readonly=False)
    uvw     = mstable.getcol('UVW')
    antable = table( os.path.join(msname, 'ANTENNA'), ack=False, readonly=True)
    baselines = int(antable.nrows() * (antable.nrows() - 1)/2. + antable.nrows())
    antable.close()
    blocksize = baselines * sbpertime
    
    if uvw.shape[0]%blocksize != 0:
        # This should match: int(uvw.shape[0]/blocksize) = xsttimes.size
        raise ValueError("\t=== Sub-bands per time  parameter 'sbertime' might be wrong ===")
    
    zenithuvw = uvw[0: blocksize]
    for it in np.arange(1, int(uvw.shape[0]/blocksize)): 
        uvw[ it*blocksize : (it+1)*blocksize ] = zenithuvw
    mstable.putcol('UVW', uvw)
    mstable.flush()
    mstable.close(); del uvw
    
    updateHist(message='UVW table re-computed towards local zenith')
    return 


# =================================================================================== #
# ----------------------------------- rephaseData ----------------------------------- #
# =================================================================================== #
def rephaseData():
    """
    """
    return


# =================================================================================== #
# ----------------------------------- updateHist ------------------------------------ #
# =================================================================================== #
def updateHist(message):
    """
    """
    return

# log = inspect.stack()
        
#         mstable = table( os.path.join(self.msfile, 'HISTORY'), ack=False, readonly=True)
#         #cli_table = mstable.getcol('CLI_COMMAND')
#         app_table = mstable.getcol('APPLICATION')
#         mes_table = mstable.getcol('MESSAGE')
#         ori_table = mstable.getcol('ORIGIN')
#         pri_table = mstable.getcol('PRIORITY')
#         tim_table = mstable.getcol('TIME')
#         mstable.flush()
#         mstable.close()

#         if app_table is None:
#             # Initialization
#             #cli_table = np.array([ log[2][4] ]) 
#             app_table = log[0][1]
#             mes_table = m
#             ori_table = log[1][3]
#             pri_table = 'INFO'
#             tim_table = Time.now().mjd * 3600. * 24.
#         else:
#             #cli_table.append( np.array([log[2][4]])  )   # {'shape': [1, 1], 'array': ['        self._makeEmptyMS()\n']}
#             app_table.append( log[0][1]  )
#             mes_table.append( m )
#             ori_table.append( log[1][3] )
#             pri_table.append( 'INFO' )
#             tim_table = np.append( tim_table, Time.now().mjd * 3600. * 24. )

#         mstable = table( os.path.join(self.msfile, 'HISTORY'), ack=False, readonly=False)
#         mstable.addrows( nrows=1 )
#         mstable.flush()
#         mstable.close()

#         mstable = table( os.path.join(self.msfile, 'HISTORY'), ack=False, readonly=False)
#         #mstable.putcol('CLI_COMMAND', cli_table)
#         mstable.putcol('APPLICATION', app_table)
#         mstable.putcol('MESSAGE'    , mes_table)
#         mstable.putcol('ORIGIN'     , ori_table)
#         mstable.putcol('PRIORITY'   , pri_table)
#         mstable.putcol('TIME'       , tim_table)
#         mstable.flush()
#         mstable.close()

#         self._taskComplete(task=log[1][3])
#         return


# =================================================================================== #
# ------------------------------------- Checks -------------------------------------- #
# =================================================================================== #
def _isMS(msname):
    """ Check that the MS exists
    """
    if not os.path.isdir(msname):
        raise IOError("\t=== MS '{}' not found ===".format(msname))
    return




