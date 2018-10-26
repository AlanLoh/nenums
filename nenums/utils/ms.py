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
import shutil

from astropy.time import Time
try:
    from pyrap.tables import table
except:
    print("\n\t=== WARNING ===")

from .astro import altaz2radec, lightSpeed
from .progressbar import ProgressBar
# from astro import altaz2radec, lightSpeed
# from progressbar import ProgressBar


__author__ = 'Alan Loh'
__copyright__ = 'Copyright 2018, nenums'
__credits__ = ['Alan Loh']
__version__ = '0.0.1'
__maintainer__ = 'Alan Loh'
__email__ = 'alan.loh@obspm.fr'
__status__ = 'Production'
__all__ = [ 'antTable',
            'emptyMS',
            'addInfos',
            'addFreq',
            'addTime',
            'addDescId',
            'addData',
            'zenithUVW',
            'rephaseData',
            'addPointing',
            'cleanDir',
            'splitMS',
            'updateHist']


# =================================================================================== #
# ------------------------------------ antTable ------------------------------------- #
# =================================================================================== #
def antTable(msname, miniarrays):
    """ Create the ANTENNA table for a Measurement Set

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **miniarrays** : np.ndarray
            Mini-Arrays from a XST file (ext: 1, key: 'noMROn')
    """
    checkAnt(miniarrays=miniarrays)
    tabname = os.path.join( os.path.dirname(msname), 'ANTENNA')
    strmas  = ",".join("'MR{}'".format(i) for i in sorted( np.squeeze(miniarrays)) ) # str list of MAs
    command = 'python /cep/lofar/nenupy/makeant.py {} "[{}]"'.format(tabname, strmas)
    os.system(command)
    return


# =================================================================================== #
# ------------------------------------- emptyMS ------------------------------------- #
# =================================================================================== #
def emptyMS(msname, start, dt, bwidth, xstsbbands):
    """ Create an empty *Measurement Set* designed for NenuFAR XST observations.
        This function calls an external program `makems` which is a LOFAR soft.

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **start** : (str, Time)
            Beginning of the observation. XST cross-correlations are computed at the local zenith.
            It helps initializing the tracked RA / Dec.
        * **dt** : float
            Time step of the XST observation, in seconds
        * **bwidth** : float
            Spectral bandwidth in Hz
        * **sbpertime** : int
            Number of sub-bands per unit of time
        * **xstsbbands** : np.ndarray
            Sub-bands of a XST file (ext: 7, key: 'xstsubband')
    """
    if os.path.isdir(msname):
        warnings.warn('Measurement Set {} already exist'.format(msname))
        return
    outpath = os.path.dirname(os.path.abspath(msname))
    if not os.path.isdir(outpath):
        os.makedirs(outpath)

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
    parset_dict['NBands']             = {'id': 0, 'val': str(xstsbbands.shape[1])}
    parset_dict['WriteAutoCorr']      = {'id': 0, 'val': 'T'}
    parset_dict['StartFreq']          = {'id': 0, 'val': '[' + str(xstsbbands[0, 0]*bwidth + bwidth/2.) + ']'}
    parset_dict['StepFreq']           = {'id': 0, 'val': str(bwidth)}
    parset_dict['NFrequencies']       = {'id': 0, 'val': str(np.unique(xstsbbands).size)}
    parset_dict['StartTime']          = {'id': 0, 'val': start.isot.replace('T', '/')}
    parset_dict['StepTime']           = {'id': 0, 'val': str(dt)}
    parset_dict['NTimes']             = {'id': 0, 'val': str(xstsbbands.shape[0])}
    parset_dict['VDSPath']            = {'id': 0, 'val': outpath}
    parset_dict['WriteImagerColumns'] = {'id': 0, 'val': "T"}
    
    # ------ Create a parset from the dictionnary ------# 
    parset = open( os.path.join(outpath, 'makems.cfg'), 'w')
    for key in parset_dict.keys():
        parset.write('{} = {}\n'.format(key, parset_dict[key]['val']))
    parset.close()
    
    # ------ Use makems to produce an empty MS ------# 
    os.system('makems {}'.format(os.path.join(outpath, 'makems.cfg')))
        
    updateHist(msname=msname, message='Creation of the Measurement Set')

    return


# =================================================================================== #
# ------------------------------------ addInfos ------------------------------------- #
# =================================================================================== #
def addInfos(msname, xstheader):
    """ Add observation standard information to the Measurement Set.
        These infos are retrieved from the XST header.

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **xstheader** : fits.header
            XST file header (ext: 0)
    """
    isMS(msname)

    mstable = table( os.path.join(msname, 'OBSERVATION'), ack=False, readonly=False)
    obse    = mstable.getcol('OBSERVER')
    proj    = mstable.getcol('PROJECT')
    sche    = mstable.getcol('SCHEDULE_TYPE')
    tele    = mstable.getcol('TELESCOPE_NAME')

    obse[0] = xstheader['CONTACT']
    proj[0] = xstheader['OBJECT']
    sche[0] = xstheader['INSTRUME']
    tele[0] = xstheader['INSTRUME']

    mstable.putcol('OBSERVER', obse)
    mstable.putcol('PROJECT', proj)
    mstable.putcol('SCHEDULE_TYPE', sche)
    mstable.putcol('TELESCOPE_NAME', tele)

    mstable.flush()
    mstable.close()

    updateHist(msname=msname, message='Observation info added, from XST fits (version {})\
        created by software version {}.'.format(xstheader['VERSION'], xstheader['SOFTWARE']))
    
    return


# =================================================================================== #
# ------------------------------------- addFreq ------------------------------------- #
# =================================================================================== #
def addFreq(msname, xstsbbands):
    """ Fill the frequency / spectral windows columns of the MS

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **xstsbbands** : np.ndarray
            Sub-bands of a XST file (ext: 7, key: 'xstsubband')
    """
    isMS(msname)

    sbpertime = xstsbbands.shape[1]
    subbands  = np.unique(xstsbbands)
    sbsize    = subbands.size
    spwtable  = table( os.path.join(msname, 'SPECTRAL_WINDOW'), ack=False, readonly=True)
    bwidth    = spwtable.getcol('EFFECTIVE_BW')[0][0]
    spwtable.close()

    # ------ Expanding SPECTRAL_WINDOW subtable ------ #
    mstable     = table( os.path.join(msname, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
    mstable.addrows( nrows=sbsize - sbpertime ) # adding rows
    mstable.flush()
    mstable.close()
    # ------ NAME column ------ #
    mstable     = table( os.path.join(msname, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
    spw_name    = mstable.getcol('NAME')
    spw_name[:] = map(str, subbands)
    mstable.putcol('NAME', spw_name) # writing "REAL" number of NenuFAR SB.
    mstable.flush()
    mstable.close()
    # ------ CHAN_FREQ column ------ #
    mstable     = table( os.path.join(msname, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
    chan_freq   = np.array([np.array([subbands[i]*bwidth + bwidth/2.]) for i in np.arange(sbsize)])
    mstable.putcol('CHAN_FREQ', chan_freq)
    mstable.flush()
    mstable.close()
    # ------ REF_FREQUENCY column ------ #
    mstable     = table( os.path.join(msname, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
    ref_freq    = chan_freq
    mstable.putcol('REF_FREQUENCY', ref_freq)
    mstable.flush()
    mstable.close()

    # ------ Taking car of other unused columns ------ #
    mstable     = table( os.path.join(msname, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
    bw_table    = np.array([ np.array([bwidth]) for i in np.arange(sbsize) ]) #  np.repeat( np.array([self.bandwidth]), self.subbands.size)
    dummy0      = np.zeros(sbsize)
    dummy1      = np.ones(sbsize)
    mstable.putcol('CHAN_WIDTH',      bw_table)
    mstable.putcol('EFFECTIVE_BW',    bw_table)
    mstable.putcol('RESOLUTION',      bw_table) 
    mstable.putcol('TOTAL_BANDWIDTH', bw_table)
    mstable.putcol('FLAG_ROW',        dummy0)
    mstable.putcol('FREQ_GROUP',      dummy0)
    mstable.putcol('IF_CONV_CHAIN',   dummy0)
    mstable.putcol('NET_SIDEBAND',    dummy0)
    mstable.putcol('NUM_CHAN',        dummy1)
    mstable.putcol('MEAS_FREQ_REF',   dummy1*5)
    mstable.flush()
    mstable.close()
    del spw_name, chan_freq, ref_freq, dummy0, dummy1, bw_table
    del mstable

    # ------ Expanding DATA_DESCRIPTION table ------ #
    mstable    = table( os.path.join(msname, 'DATA_DESCRIPTION'), ack=False, readonly=False)
    mstable.addrows( nrows=sbsize - sbpertime ) # adding rows
    # ------ Inserting index to map DATA_DESC_ID index to the SPECTRAL_WINDOW index ------ #
    spwid_table    = mstable.getcol('SPECTRAL_WINDOW_ID')
    spwid_table[:] = np.arange(spwid_table.size)
    mstable.putcol('SPECTRAL_WINDOW_ID', spwid_table)
    dummy          = np.zeros(spwid_table.size) # create dummy array
    # ------ Filling all columns ------- #
    mstable.putcol('FLAG_ROW',           dummy)
    mstable.putcol('POLARIZATION_ID',    dummy)
    mstable.putcol('SPECTRAL_WINDOW_ID', spwid_table)
    mstable.flush()
    mstable.close()
    del mstable, spwid_table, dummy
    
    updateHist(msname=msname, message='SPECTRAL_WINDOW and DATA_DESCRIPTION tables written.')
    return


# =================================================================================== #
# ------------------------------------- addTime ------------------------------------- #
# =================================================================================== #
def addTime(msname, xsttime):
    """ Fill the TIME column of the MS based on XST times (in Julian Days)

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **xsttime** : np.ndarray
            Time array from a XST file (in JD)
    """
    isMS(msname)

    antable   = table( os.path.join(msname, 'ANTENNA'), ack=False, readonly=True)
    baselines = int(antable.nrows() * (antable.nrows() - 1)/2. + antable.nrows())
    antable.close()

    mstable    = table(msname, ack=False, readonly=False)
    mstime     = mstable.getcol('TIME')
    umstime    = np.unique(mstime)
    dt         = umstime[1] - umstime[0]
    time_table = Time(xsttime, format='jd').mjd * 3600. * 24. + (dt/2.) # MJD in sec (MS standards)
    sbpertime  = int(mstime.size / baselines / time_table.size)
    time_table = np.repeat(time_table, baselines * sbpertime)
    mstable.putcol('TIME', time_table) 
    mstable.putcol('TIME_CENTROID', time_table)  
    mstable.flush()
    mstable.close()

    updateHist(msname=msname, message='TIME table filled.')
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
    isMS(msname)

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

    updateHist(msname=msname, message='DATA_DESC_ID table filled.')
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
    isMS(msname)

    antable   = table( os.path.join(msname, 'ANTENNA'), ack=False, readonly=True)
    baselines = int(antable.nrows() * (antable.nrows() - 1)/2. + antable.nrows())
    antnames  = antable.getcol('NAME')
    usedmas   = np.array([ int(mr.replace('MR', '')) for mr in antnames ])
    antable.close()

    mstable    = table(msname, ack=False, readonly=False)
    data_ms    = mstable.getcol('DATA')
    ant1_table = mstable.getcol('ANTENNA1')
    ant2_table = mstable.getcol('ANTENNA2')

    bar = ProgressBar(valmax=baselines, title='Filling MS data from the XST cross-correlations')
    for ant1 in builtma:
        for ant2 in builtma[ant1:]:
            if (ant1 not in usedmas) or (ant2 not in usedmas):
                continue        
            baseline = np.sort( [ant1, ant2] )[::-1] # descending order
            # These formulas were defined to fit the table of detail_FITS_00_05.pdf (kinda triangular series)
            MA1X_MA2X = baseline[0]*2    *(baseline[0]*2+1)/2   + 2*baseline[1]
            MA1X_MA2Y = baseline[0]*2    *(baseline[0]*2+1)/2+1 + 2*baseline[1] # it will be the same as next one for auto-corr
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
            data_ms[clever_index, :, :] = np.reshape(data, (data.shape[0]*data.shape[1], 1, data.shape[2])) # (ntimes*nsubpertime, nchan, npol)
            bar.update()

    mstable.putcol('DATA', data_ms) 
    mstable.flush()
    mstable.close()

    checkSum(msname=msname)

    updateHist(msname=msname,  message='DATA table filled with XST cross-correlations.')
    return


# =================================================================================== #
# ----------------------------------- zenithUVW ------------------------------------- #
# =================================================================================== #
def zenithUVW(msname):
    """ Re-compute the UVW table of the Measurement Set.
        `emptyMS()` automatically creates a UVW table as if initial RA Dec (zenith at start) was tracked.
        However, XST are cross-correlations fixed at the local zenith.
        Therefore, the correct UVW table must be copied from the first time step.

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
    """
    isMS(msname)
    
    mstable = table(msname, ack=False, readonly=False)
    uvw     = mstable.getcol('UVW')
    antable = table( os.path.join(msname, 'ANTENNA'), ack=False, readonly=True)
    baselines = int(antable.nrows() * (antable.nrows() - 1)/2. + antable.nrows())
    antable.close()
    
    mstable    = table(msname, ack=False, readonly=True)
    mstime     = mstable.getcol('TIME')
    sbpertime  = int(mstime.size / baselines / np.unique(mstime).size)

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
    
    updateHist(msname=msname, message='UVW table re-computed towards local zenith')
    return 


# =================================================================================== #
# ----------------------------------- rephaseData ----------------------------------- #
# =================================================================================== #
def rephaseData(msname, xsttime, ra_center, dec_center):
    """ XST Cross-Correlations are phased at the local zenith.
        UVW are supposed to also be computed at the local zenith (`zenithUVW()`).
        In order to produce a tracking observation, this function recompute UVW
        and re-phase the visibilities towards a new phase center (`ra_center`, `dec_center`)

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **xsttime** : np.ndarray
            Time array from a XST file (in JD)
        * **ra_center** : float
            Right-Ascension of the tracked target (in radians)
        * **dec_center** : float
            Declination of the tracked target (in radians)
    """
    isMS(msname)

    # ------ Original phase centers ------ #
    az_list  = np.repeat(0. , xsttime.size) # degrees
    alt_list = np.repeat(90., xsttime.size)
    ra_zenith, dec_zenith = altaz2radec(az_list, alt_list, Time(xsttime, format='jd'))

    # ------ Get data to transform ------
    mstable   = table(os.path.join(msname), ack=False, readonly=True)
    uvw       = mstable.getcol('UVW')
    data      = mstable.getcol('DATA')
    spwdataid = mstable.getcol('DATA_DESC_ID')  # Nvis*Ntimes*Kspw
    spwtable  = table(os.path.join(msname, 'SPECTRAL_WINDOW'), ack=False, readonly=True)
    chanfreq  = spwtable.getcol('CHAN_FREQ')   # Nvis*Ntimes*Kspw x Nchan x 4   in Hz
    desctable = table( os.path.join(msname, 'DATA_DESCRIPTION'), ack=False, readonly=True)
    spwid     = desctable.getcol('SPECTRAL_WINDOW_ID')  # Total list of spectral window ID
    mstable.close()
    spwtable.close()
    desctable.close()

    # ------ Properties ------ #
    antable   = table( os.path.join(msname, 'ANTENNA'), ack=False, readonly=True)
    baselines = int(antable.nrows() * (antable.nrows() - 1)/2. + antable.nrows())
    antable.close()

    mstable    = table(msname, ack=False, readonly=True)
    mstime     = mstable.getcol('TIME')
    sbpertime  = int(mstime.size / baselines / np.unique(mstime).size)
    mstable.close()

    # ------ Wavelength ------
    frequency  = np.take(chanfreq, spwid[spwdataid], axis=0) # Nvis*Ntimes*Nspw*Nchan extract of list of channel frequencies for the relevant SPWID
    wavelength = lightSpeed() / frequency

    # ------ Re-phase ------
    # Going from zenith ra / dec (time dependent) towards new ra_center, dec_center
    final_trans, w_new = rotMatrix(ra_center, dec_center)
    blocksize = baselines * sbpertime
    bar = ProgressBar(valmax=ra_zenith.size, title='Re-phasing the data towards RA='+str(np.degrees(ra_center))+' DEC='+str(np.degrees(dec_center)))
    for itime in np.arange(ra_zenith.size):
        idx1 = itime     * blocksize # start index
        idx2 = (itime+1) * blocksize
        original_trans, w_old = rotMatrix(ra_zenith[itime], dec_zenith[itime]) 
        total_trans           = np.dot(final_trans.T, original_trans)
        if itime > 1:
            # Check that original UVW tables are identical for every time step
            if not np.array_equal(tmpuvw, uvw[idx1:idx2]):
                raise ValueError("\t=== Initial UVW are not equal at every time frame ===")
        tmpuvw = uvw[idx1:idx2].copy()
        phase  = np.dot( np.dot( (w_old-w_new).T, original_trans) , uvw[idx1:idx2].T )

        # Updating UVW
        uvw[idx1:idx2] = np.dot(uvw[idx1:idx2], total_trans.T)
        nrows, nchan, npol = data[idx1:idx2].shape
        for chan in np.arange(nchan):
            dphi = np.exp( phase * 2 * np.pi * 1j / wavelength[idx1:idx2, chan])
            for pol in np.arange(npol):
                # Data rephasing
                data[idx1:idx2, chan, pol] = data[idx1:idx2, chan, pol] * dphi
        bar.update()

    # ------ Writing Output ------ #
    mstable = table(msname, ack=False, readonly=False)
    mstable.putcol('UVW', uvw)
    mstable.flush()
    mstable.close(); del uvw

    mstable = table(msname, ack=False, readonly=False)
    mstable.putcol('DATA', data)
    mstable.flush()
    mstable.close(); del data

    updateHist(msname=msname, message='UVW and DATA rephased re-computed to RA='+str(np.degrees(ra_center))+' DEC='+str(np.degrees(dec_center)))
    return


# =================================================================================== #
# ----------------------------------- addPointing ----------------------------------- #
# =================================================================================== #
def addPointing(msname, ra_center, dec_center):
    """ Add pointing informations to the MS

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
         * **ra_center** : float
            Right-Ascension of the tracked target (in radians)
        * **dec_center** : float
            Declination of the tracked target (in radians)
    """
    isMS(msname)

    # ------ Changing Phase Center in FIELD subtable ------ # 
    mstable   = table( os.path.join(msname, 'FIELD'), ack=False, readonly=False)
    delay_tab = mstable.getcol('DELAY_DIR')
    phase_tab = mstable.getcol('PHASE_DIR')
    refer_tab = mstable.getcol('REFERENCE_DIR')
    delay_tab[0, 0, 0] = ra_center # radians !
    delay_tab[0, 0, 1] = dec_center
    phase_tab[0, 0, 0] = ra_center
    phase_tab[0, 0, 1] = dec_center
    refer_tab[0, 0, 0] = ra_center  
    refer_tab[0, 0, 1] = dec_center 
    mstable.putcol('DELAY_DIR', delay_tab)
    mstable.putcol('PHASE_DIR', phase_tab)
    mstable.putcol('REFERENCE_DIR', refer_tab)
    mstable.flush()
    mstable.close()

    # ------ Changing POINTING subtable ------ # 
    mstable   = table( os.path.join(msname, 'POINTING'), ack=False, readonly=False)
    dir_tab   = mstable.getcol('DIRECTION')
    targ_tab  = mstable.getcol('TARGET')
    track_tab = mstable.getcol('TRACKING')
    dir_tab[:,  0, 0] = np.repeat(ra_center,  dir_tab[:,  0, 0].size)
    dir_tab[:,  0, 1] = np.repeat(dec_center, dir_tab[:,  0, 1].size)
    targ_tab[:, 0, 0] = np.repeat(ra_center,  targ_tab[:, 0, 0].size)
    targ_tab[:, 0, 1] = np.repeat(dec_center, targ_tab[:, 0, 1].size)
    track_tab[:] = np.ones(track_tab.size, dtype=bool)
    mstable.putcol('DIRECTION', dir_tab)
    mstable.putcol('TARGET', targ_tab)
    mstable.putcol('TRACKING', track_tab)
    mstable.flush()
    mstable.close()

    updateHist(msname=msname, message='Pointing informations updated')
    return


# =================================================================================== #
# ------------------------------------ cleanDir ------------------------------------- #
# =================================================================================== #
def cleanDir(msname):
    """ Remove unused sub-items created during the MS construction process

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
    """
    msname = os.path.abspath(msname)
    mspath = os.path.dirname(msname)
    ant_table = os.path.join(mspath, 'ANTENNA')
    try:
        shutil.rmtree(ant_table, ignore_errors=True)
        os.rename( os.path.join(mspath, 'makems.cfg'), os.path.join(mspath, os.path.basename(msname)[:-3] + '_makems.cfg'))
        os.remove( msname + '.gds' )
        os.remove( msname + '.vds' )
    except:
        print("\t=== WARNING: error trying to clean {} from unused files. ===".format(mspath))
    return


# =================================================================================== #
# ------------------------------------- splitMS ------------------------------------- #
# =================================================================================== #
def splitMS(msname, remove=False):
    """ Split the MS by sub-bands, as LOFAR

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **remove** : bool, optional
            Remove the Measurement Set once it has correctly been split (default = False)
    """
    mspath = os.path.dirname(msname)

    freqt = table( os.path.join(msname, 'SPECTRAL_WINDOW'), ack=False, readonly=True)
    names = freqt.getcol('NAME')
    freqs = freqt.getcol('REF_FREQUENCY')
    subbands = freqt.nrows()
    freqt.close()
    desct = table(os.path.join(msname, 'DATA_DESCRIPTION'), ack=False, readonly=True)
    spws  = desct.getcol('SPECTRAL_WINDOW_ID')
    desct.close()

    mstable = table(msname, ack=False, readonly=True)
    for spw, name, freq in zip(spws, names, freqs):
        output = msname.split('.ms')[0] + '_SB' + str(name) + '.ms'
        newms  = mstable.query('DATA_DESC_ID == {}'.format(spw))
        newms.copy(output, deep=True)
        newms.close()
        # Update SPECTRAL_WINDOW and DATA_DESCRIPTION tables to match the frequency selection
        ftable = table( os.path.join(output, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
        ftable.removerows(rownrs=np.arange(subbands)[np.arange(subbands) != spw])
        ftable.flush()
        ftable.close()
        mstab  = table(output, ack=False, readonly=False)
        datad  = mstab.getcol('DATA_DESC_ID')
        datad  = np.zeros( datad.shape ) # every element is labelled 'SPW=0' now...
        mstab.putcol('DATA_DESC_ID', datad)
        mstab.flush()
        mstab.close()
        dtable = table( os.path.join(output, 'DATA_DESCRIPTION'), ack=False, readonly=False)
        dtable.removerows(rownrs=np.arange(subbands)[np.arange(subbands) != spw])
        spwid = dtable.getcol('SPECTRAL_WINDOW_ID')
        spwid[0] = 0 # set the SPW index to 0 since this is the only frequency
        dtable.putcol('SPECTRAL_WINDOW_ID', spwid)
        dtable.flush()
        dtable.close()

        updateHist(msname=output, message='Split from MS {}'.fotmat(msname))

        print("\t=== Measurement set {} written. ===".format(output))

    mstable.close()

    # ------ Remove the MS ------ #
    if remove:
        sbms = glob.glob( os.path.join(mspath, '*_SB*.ms') )
        if len(sbms) != subbands:
            print("\t=== WARNING: could not find {} sub-band MSs in {}. ===".format(subbands, mspath))
        else:
            try:
                shutil.rmtree(msname, ignore_errors=True)
            except:
                print("\t=== WARNING: impossible to remove {} ===".format(msname))
        pass
    return


# =================================================================================== #
# ----------------------------------- updateHist ------------------------------------ #
# =================================================================================== #
def updateHist(msname, message):
    """ Keep track of the modifications.
        Store everything into the HISTORY table of the Measurment set

        Parameters
        ----------
        * **msname** : str
            Name of the Measurement Set
        * **message** : str
            Message to be kept in HISTORY
    """
    log = inspect.stack()
        
    mstable = table( os.path.join(msname, 'HISTORY'), ack=False, readonly=True)
    #cli_table = mstable.getcol('CLI_COMMAND')
    app_table = mstable.getcol('APPLICATION')
    mes_table = mstable.getcol('MESSAGE')
    ori_table = mstable.getcol('ORIGIN')
    pri_table = mstable.getcol('PRIORITY')
    tim_table = mstable.getcol('TIME')
    mstable.flush()
    mstable.close()

    if app_table is None:
        # Initialization
        #cli_table = np.array([ log[2][4] ]) 
        app_table = log[0][1]
        mes_table = message
        ori_table = log[1][3]
        pri_table = 'INFO'
        tim_table = Time.now().mjd * 3600. * 24.
    else:
        #cli_table.append( np.array([log[2][4]])  )   # {'shape': [1, 1], 'array': ['        self._makeEmptyMS()\n']}
        app_table.append( log[0][1]  )
        mes_table.append( message )
        ori_table.append( log[1][3] )
        pri_table.append( 'INFO' )
        tim_table = np.append( tim_table, Time.now().mjd * 3600. * 24. )

    mstable = table( os.path.join(msname, 'HISTORY'), ack=False, readonly=False)
    mstable.addrows( nrows=1 )
    mstable.flush()
    mstable.close()

    mstable = table( os.path.join(msname, 'HISTORY'), ack=False, readonly=False)
    #mstable.putcol('CLI_COMMAND', cli_table)
    mstable.putcol('APPLICATION', app_table)
    mstable.putcol('MESSAGE'    , mes_table)
    mstable.putcol('ORIGIN'     , ori_table)
    mstable.putcol('PRIORITY'   , pri_table)
    mstable.putcol('TIME'       , tim_table)
    mstable.flush()
    mstable.close()

    # self._taskComplete(task=log[1][3])
    return


# =================================================================================== #
# ----------------------------------- rotMatrix ------------------------------------- #
# =================================================================================== #
def rotMatrix(ra, dec):
    """ Compute a 3D rotation matrix
        
        Parameters
        ----------
        * **ra** : float
            Right Ascension (in radians)
        * **dec** : float
            Declination (in radians)

        Returns
        -------
        * **rot** : np.ndarray
            Rotation matrix
        * **w** : np.ndarray
            W vector
    """
    w = np.array([[  np.sin(ra)*np.cos(dec) ,  np.cos(ra)*np.cos(dec) , np.sin(dec) ]]).T
    v = np.array([[ -np.sin(ra)*np.sin(dec) , -np.cos(ra)*np.sin(dec) , np.cos(dec) ]]).T
    u = np.array([[  np.cos(ra)             , -np.sin(ra)             , 0.          ]]).T
    rot_matrix = np.concatenate([u, v, w], axis=-1)
    return rot_matrix, w


# =================================================================================== #
# ------------------------------------- Checks -------------------------------------- #
# =================================================================================== #
def isMS(msname):
    """ Check that the MS exists
    """
    if not os.path.isdir(msname):
        raise IOError("\t=== MS '{}' not found ===".format(msname))
    return


def checkAnt(miniarrays):
    """ Check that all required Mini-Arrays are in the location file
        
        Parameters
        ----------
        * **miniarrays** : np.ndarray
            Mini-Arrays from a XST file (ext: 1, key: 'noMROn')
    """
    modulepath = os.path.dirname( os.path.realpath(__file__) )
    rfile    = open(os.path.join(modulepath, 'locationsSB2.dat'), 'r')
    mrinfile = []
    for line in rfile:
        if "DBNAME='NENUFAR'" in line:
            mrinfile.append( int(line.split("'")[1][2:]) )
    for ma in miniarrays:
        if ma not in mrinfile:
            raise ValueError("\t=== MA {} not in {} ===".format(ma, os.path.join(modulepath, 'locationsSB2.dat')))
    return


def checkSum(msname):
    """ Check that the sum of all auto-correlations imag( XY ) + imag( YX ) is 0.
        This method is completely independent from the way data are stored in the MS.
    """
    mstable = table(msname, ack=False, readonly=True)
    mstable = mstable.query('ANTENNA1 == ANTENNA2') # get only auto-correlations
    data   = mstable.getcol('DATA')
    mstable.close()
    auto_sum = np.sum(data[:, : , 1].imag + data[:, :, 2 ].imag)
    if auto_sum != 0.:
        raise ValueError("\t=== WARNING: Something's wrong... Failed to verify _checkSum() test! ===")
    else:
        # Everything's alright ;-)
        pass
    return

