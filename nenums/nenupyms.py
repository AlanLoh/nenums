#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = "Alan Loh, Julien Girard"
__credits__ = ["Alan Loh", "Julien Girard"]
__version__ = "0.0.1"
__maintainer__ = "Alan Loh"
__email__ = "alan.loh@obspm.fr"
__status__ = "Production"

# --------------------------------------------------------------------------------------------- #
# ---------------------------------------- TO DO LIST  ---------------------------------------- #
# - fonction de clean du repertoire -- DONE
# - re-calcul des UVW -- DONE
# - re-phasage des donnees apres calcul correct uvw  -- DONE
# - test de somme des autoccorelations im + re pour etre sur d avoir 0 -- DONE
# - definition automatique des spw en fonction des sous bandes visitees
# - faire une table history du ms -- DONE
# - ajouter info pointing ou tracking dans les headers des ms
# - rajouter une progress bar sur la boucle baseline dans data -- DONE
# - nom du MS en fonction de celui du XST type 'NENUFAR_20180601_021250.ms' -- DONE
# - retirer les tables antennes et mettre les logs dans un dossier
# - ajouter un truc pour supprimer le ms après un split (le laisser en option peut etre)
# - mettre une option pour copier DATA dans CORRECTED_DATA -- DONE
# - compresser le tout en .tar.gz
# --------------------------------------------------------------------------------------------- #


import nenupy
import os
import shutil
import tarfile
import inspect
import numpy as np
from pyrap.tables import table, addImagingColumns
from astropy.time import Time
from astropy import units as u
from astropy import coordinates as coord
from astropy import constants as const
from progress.bar import ChargingBar

class NenuFits(nenupy.Read):
    def __init__(self, obsfile, savepath='', **kwargs):
        nenupy.Read.__init__(self, obsfile, **kwargs)
        if self.type != 'XST':
            raise IOError("\t=== Observation file needs to be a valid XST fits. ===")
        self._antlocation = '/cep/lofar/nenupy/locationsSB2.dat'
        #self.nancay       = coord.EarthLocation(lat=47.3759883*u.deg, lon=2.1927188*u.deg)
        self.nancay       = coord.EarthLocation(lat=47.376511*u.deg, lon=2.1924002*u.deg)
        self.savepath     = savepath
        self.subpertime   = 16 # number os sub-bands per unit of time
        self.nb_baseline  = self.miniarrays.size * (self.miniarrays.size - 1)/2 + self.miniarrays.size
        self.subbands     = np.unique(self.xstsubb)
        #self.spec_windows = np.arange(self.subbands.size) # SPW indices
        self.do_tracking  = True # Re-phase
        self.cp_datacol   = False # Coypy the data column to the CORRECTED_DATA
        self.forceradec   = None # otherwise (ra, dec) in deg
    
    # --------------------------------------------------------------------------------------------- #
    # ----------------------------------- Attributes Checking  ------------------------------------ #
    @property
    def savepath(self):
        """ Path to save a Measurement Set
        """
        return self._savepath
    @savepath.setter
    def savepath(self, s):
        s = os.path.abspath(s)
        if os.path.isdir(s):
            pass
        else:
            try:
                os.mkdir(s)
                
            except:
                raise IOError("\t=== Impossible to create directory {} ===".format(s))
        self._savepath = s
        self.msfile    = os.path.join(self._savepath, 'NENUFAR_XST_{}.ms'.format(self.obsname))
        return
    # ------ Specific getters ------ #
    @property
    def spec_windows(self):
        """ Spectral windows definitions
        """
        # ------ Detect if contiguous sub-bands ------ #
        groups = np.split(self.subbands, np.where(np.diff(self.subbands) != 1)[0] + 1)
        if np.unique( [len(groups[i]) for i in xrange(len(groups))] ).size == 1:
            # The groups contain the same amount of sub-bands
            spw = np.arange( len(groups) )
        else:
            # Create as many SPW as sub-bands
            spw = np.arange(self.subbands.size)
        return spw

    # --------------------------------------------------------------------------------------------- #
    # ----------------------------------------- Methods ------------------------------------------- #
    def createMs(self, sbsplit=False):
        """ Create a Measurement Set
            import nenupyms; a=nenupyms.NenuFits('20180601_021250_XST.fits'); a.createMs()
            import nenupyms; a=nenupyms.NenuFits('/databf/nenufar/20180615_135300_20180615_135910_ZENITH_TEST_MROFF/20180615_135300_XST.fits'); a.savepath='/data/loh/NenuFAR/TEST_MS/20180615_MR_OFF'; a.createMs()
            import nenupyms; a=nenupyms.NenuFits('/data/loh/NenuFAR/Test_MS_labo/20180605_125000_XST.fits', '/data/loh/NenuFAR/Test_MS_labo/FINAL_MS'); a.createMs()
            import nenupyms; a=nenupyms.NenuFits('/databf/nenufar/20180602_020900_20180602_040920_CYGA_TRACKING_CALIBI/20180602_020900_XST.fits', '/data/loh/NenuFAR/TEST_MS/20180602_CygA'); a.createMs()
            import nenupyms; a=nenupyms.NenuFits('/databf/nenufar/20180603_023500_20180603_033520_CYGA_TRACKING_CALIBI/20180603_023500_XST.fits', '/data/loh/NenuFAR/TEST_MS/20180603_CygA_TEST'); a.createMs()
            import nenupyms; a=nenupyms.NenuFits('/databf/nenufar/20180603_055830_20180603_065850_CASA_TRACKING_CALIBI/20180603_055830_XST.fits', '/data/loh/NenuFAR/A_TEAM_IMAGING/20180603_CasA_CYGA'); a.forceradec = (299.86815263, 40.73391583); a.createMs()
        """
        self._makeAntTab()
        self._makeEmptyMS()
        self._fillObsinfos()
        self._fillFreq()
        self._fillData()
        self._fillUVW()
        self._fillPointing()
        self._addCorrected()
        self._cleanDirectory()
        if sbsplit:
            self._splitMS()
        return

    # --------------------------------------------------------------------------------------------- #
    # --------------------------------------- MS Creation ----------------------------------------- #
    def _makeAntTab(self):
        """ Construction of the ANTENNA table
        """
        self._checkAnt()
        tabname = os.path.join( os.path.dirname( os.path.abspath(self.msfile) ), 'ANTENNA')
        strmas  = ",".join("'MR{}'".format(i) for i in sorted(self.miniarrays) ) # str list of MAs
        command = 'python /cep/lofar/nenupy/makeant.py {} "[{}]"'.format(tabname, strmas)
        os.system(command)
        return
    def _makeEmptyMS(self):
        """ Create an empty Measurement Set
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
    def _fillObsinfos(self):
        """ Fill the OBSERVATION subtable
        """
        mstable   = table( os.path.join(self.msfile, 'OBSERVATION'), ack=False, readonly=False)
        observer  = mstable.getcol('OBSERVER')
        project   = mstable.getcol('PROJECT')
        schedule  = mstable.getcol('SCHEDULE_TYPE')
        telescope = mstable.getcol('TELESCOPE_NAME')

        observer[0]  = self.header['CONTACT']
        project[0]   = self.header['OBJECT']
        schedule[0]  = self.header['INSTRUME']
        telescope[0] = self.header['INSTRUME']

        mstable.putcol('OBSERVER', observer)
        mstable.putcol('PROJECT', project)
        mstable.putcol('SCHEDULE_TYPE', schedule)
        mstable.putcol('TELESCOPE_NAME', telescope)

        mstable.flush()
        mstable.close()

        del observer, project, schedule, telescope
        del mstable
        self._updateHistory(m='Observation info added, from XST fits (version {}) created by software version {}.'.format(self.header['VERSION'], self.header['SOFTWARE']))
        return
    def _fillFreq(self):
        """ Addd the frequency column to the MS
        """
        # ------ Expanding SPECTRAL_WINDOW subtable ------ #
        mstable     = table( os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
        mstable.addrows( nrows=self.subbands.size - self.subpertime ) # adding rows
        mstable.flush()
        mstable.close()
        # ------ NAME column ------ #
        mstable     = table( os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
        spw_name    = mstable.getcol('NAME')
        spw_name[:] = map(str, self.subbands)
        mstable.putcol('NAME', spw_name) # writing "REAL" number of NenuFAR SB.
        mstable.flush()
        mstable.close()
        # ------ CHAN_FREQ column ------ #
        mstable     = table( os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
        chan_freq   = np.array([np.array([self.subbands[i] * self.bandwidth + self.bandwidth/2.]) for i in xrange(self.subbands.size)])
        mstable.putcol('CHAN_FREQ', chan_freq)
        mstable.flush()
        mstable.close()
        # ------ REF_FREQUENCY column ------ #
        mstable     = table( os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
        ref_freq    = chan_freq
        mstable.putcol('REF_FREQUENCY', ref_freq)
        mstable.flush()
        mstable.close()

        # ------ Taking car of other unused columns ------ #
        mstable     = table( os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
        bw_table    = np.array([ np.array([self.bandwidth]) for i in xrange(self.subbands.size) ]) #  np.repeat( np.array([self.bandwidth]), self.subbands.size)
        dummy0      = np.zeros(self.subbands.size)
        dummy1      = np.ones(self.subbands.size)
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
        mstable    = table( os.path.join(self.msfile, 'DATA_DESCRIPTION'), ack=False, readonly=False)
        mstable.addrows( nrows=self.subbands.size - self.subpertime ) # adding rows
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
        self._updateHistory(m='SPECTRAL_WINDOW and DATA_DESCRIPTION tables written.')
        return
    def _fillData(self):
        """ 
        """
        # ------ Time ------ #
        mstable    = table(self.msfile, ack=False, readonly=False)
        time_table = Time(self.xsttime, format='jd').mjd * 3600. * 24. + (self.steptime/2.) # MJD in sec (MS standards), centered via steptime
        time_table = np.repeat(time_table, self.nb_baseline * self.subpertime)
        mstable.putcol('TIME', time_table)  
        mstable.flush()
        mstable.close()

        # ------ Data Desc Id ------ #
        mstable = table(self.msfile, ack=False, readonly=False)
        desc_table   = mstable.getcol('DATA_DESC_ID') # nvis * ntime
        spw_idx_time = np.searchsorted(self.subbands, self.xstsubb.ravel()) # associate the SPW index to observed Sub-Band
        # Make sure that unused SPW are not indexed
        mask         = spw_idx_time < self.subbands.size
        mask[mask]   = (self.subbands[spw_idx_time[mask]] == self.xstsubb.ravel()[mask])
        spw_idx_time = spw_idx_time[mask] # SPW index per unit of time
        spw_idx_time = np.repeat(spw_idx_time, self.nb_baseline) # repeat for the number of visibilities

        mstable.putcol('DATA_DESC_ID', spw_idx_time)  
        mstable.flush()
        mstable.close()
        del mstable, spw_idx_time, mask

        # ------ DATA !! ------ #
        mstable = table(self.msfile, ack=False, readonly=False)
        data_ms = mstable.getcol('DATA')
        ant1_table = mstable.getcol('ANTENNA1')
        ant2_table = mstable.getcol('ANTENNA2')

        bar = ChargingBar('Filling MS data from the XST cross-correlations', max=self.nb_baseline)
        for ant1 in self._allma:
            for ant2 in self._allma[ant1:]:
                if (ant1 not in self.miniarrays) or (ant2 not in self.miniarrays):
                    continue
                else:
                    # Get the time x frequency x correl for the current baseline
                    baseline = np.sort( [ant1, ant2] )[::-1] # descending order
                    # ------ Indices ------
                    # These formulas were defined to fit the table of detail_FITS_00_05.pdf (it a kind of triangular series)
                    MA1X_MA2X = baseline[0]*2*(baseline[0]*2+1)/2       + 2*baseline[1]
                    MA1X_MA2Y = baseline[0]*2*(baseline[0]*2+1)/2+1     + 2*baseline[1] # it will be the same as next one for auto-corr
                    MA1Y_MA2X = (baseline[0]*2+1)*(baseline[0]*2+2)/2   + 2*baseline[1]
                    MA1Y_MA2Y = (baseline[0]*2+1)*(baseline[0]*2+2)/2+1 + 2*baseline[1]
                    index     = np.array([MA1X_MA2X, MA1X_MA2Y, MA1Y_MA2X, MA1Y_MA2Y])
                    data      = self.xstdata[:, :, index] # (time, frequency, 4)
                    if ant1 == ant2:
                        # Index 1 (XY) is wrong for auto-correlations with the above algorithm,
                        # doesn't matter, we replace the data by the conjugate of YX which exists anyway.
                        data[:, :, 1] = np.conjugate(data[:, :, 2]) 
                    ant1_index   = np.where(self.miniarrays == ant1)[0][0]
                    ant2_index   = np.where(self.miniarrays == ant2)[0][0]
                    clever_index = (ant1_table == ant1_index) & (ant2_table == ant2_index)
                    data_ms[clever_index, :, :] = np.reshape(data, (data.shape[0]*data.shape[1], 1, data.shape[2]))
                    bar.next()
        bar.finish()

        # ------- Write the data ------ # 
        mstable.putcol('DATA', data_ms) 
        mstable.flush()
        mstable.close()

        # ------- Check that the data are good to go ------ #
        self._checkSum()

        self._updateHistory(m='DATA table filled with XST cross-correlations.')
        return
    def _fillUVW(self):
        """ Re-compute UVW
            Data are zenith-phased
            The first uvw time step is good (makems computes it on the RA DEC assocoiated with the zenith at the start of obs)
        """
        #uvwfile = os.path.join(self.savepath, 'uvw.npy') # file to save the first time step, in case of error
        mstable = table(self.msfile, ack=False, readonly=False)
        uvw = mstable.getcol('UVW')
        blocksize = self.nb_baseline * self.subpertime
        for itime in xrange(1, self.xsttime.size):
            idx1 = itime     * blocksize
            idx2 = (itime+1) * blocksize
            uvw[idx1:idx2] = uvw[0: blocksize]
        # uvw[:] = np.tile(uvw[0: blocksize], (self.xsttime.size, 1) ) # way cooler but slower :-( #WeirdPython
        mstable.putcol('UVW', uvw)
        mstable.flush()
        mstable.close(); del uvw
        self._updateHistory(m='UVW table re-computed towards local zenith')
        return 
    def _fillPointing(self):
        """ Re-compute UVW to account for tracking 
            Re-phase visibility towards a new phase center (ra, dec)
        """
        if not self.do_tracking:
            return
        def rotMatrix(r, d):
            """ r: ra in radians
                d: dec in radians
            """
            w = np.array([[  np.sin(r)*np.cos(d) ,  np.cos(r)*np.cos(d) , np.sin(d) ]]).T
            v = np.array([[ -np.sin(r)*np.sin(d) , -np.cos(r)*np.sin(d) , np.cos(d) ]]).T
            u = np.array([[  np.cos(r)           , -np.sin(r)           , 0.        ]]).T
            rot_matrix = np.concatenate([u, v, w], axis=-1)
            return rot_matrix, w
        
        # ------ Phase center positions ------
        ra_z, dec_z = self._radecTime() # Original RA DEC phase center (radians), same dimension as self.xsttime
        if self.forceradec is None:
            ra_f, dec_f = self._azel2radec(self._azlistdig[0][0], self._ellistdig[0][0], Time(self.xsttime[0], format='jd')) # Final phase center
        else:
            try:
                ra_f, dec_f = np.radians(self.forceradec)
            except:
                raise ValueError("\t=== Error, self.forceradec must be (ra, dec) in degrees, found '({})'. ===".format(self.forceradec))

        # ------ Get data to transform ------
        mstable   = table(os.path.join(self.msfile), ack=False, readonly=True)
        uvw       = mstable.getcol('UVW')
        data      = mstable.getcol('DATA')
        spwdataid = mstable.getcol('DATA_DESC_ID')  # Nvis*Ntimes*Kspw
        spwtable  = table(os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=True)
        chanfreq  = spwtable.getcol('CHAN_FREQ')   # Nvis*Ntimes*Kspw x Nchan x 4   in Hz
        desctable = table( os.path.join(self.msfile, 'DATA_DESCRIPTION'), ack=False, readonly=True)
        spwid     = desctable.getcol('SPECTRAL_WINDOW_ID')  # Total list of spectral window ID
        mstable.close()
        spwtable.close()
        desctable.close()

        # ------ Wavelength ------
        frequency  = np.take(chanfreq, spwid[spwdataid], axis=0) # Nvis*Ntimes*Nspw*Nchan extract of list of channel frequencies for the relevant SPWID
        wavelength = const.c.value / frequency

        # ------ Re-phase ------
        # Going from zenith ra / dec (time dependent) towards new ra_f, dec_f
        final_trans, w_new = rotMatrix(ra_f, dec_f)
        blocksize = self.nb_baseline * self.subpertime
        bar = ChargingBar('Re-phasing the data towards RA='+str(np.degrees(ra_f))+' DEC='+str(np.degrees(dec_f)), max=ra_z.size)
        for itime in xrange(ra_z.size):
            idx1 = itime     * blocksize # start index
            idx2 = (itime+1) * blocksize
            original_trans, w_old = rotMatrix(ra_z[itime], dec_z[itime]) 
            total_trans           = np.dot(final_trans.T, original_trans)
            if itime > 1:
                # Check that original UVW tables are identical for every time step
                if not np.array_equal(tmpuvw, uvw[idx1:idx2]):
                    raise ValueError("\t=== Initial UVW are not equal at every time frame ===")
            tmpuvw = uvw[idx1:idx2].copy()
            phase                 = np.dot( np.dot( (w_old-w_new).T, original_trans) , uvw[idx1:idx2].T )
            # Updating UVW
            uvw[idx1:idx2]        = np.dot(uvw[idx1:idx2], total_trans.T)
            nrows, nchan, npol = data[idx1:idx2].shape
            for chan in xrange(nchan):
                dphi = np.exp( phase * 2 * np.pi * 1j / wavelength[idx1:idx2, chan])
                for pol in xrange(npol):
                    # Data rephasing
                    data[idx1:idx2, chan, pol] = data[idx1:idx2, chan, pol] * dphi
            bar.next()
        bar.finish()

        # ------ Writing Output ------ #
        mstable = table(self.msfile, ack=False, readonly=False)
        mstable.putcol('UVW', uvw)
        mstable.flush()
        mstable.close(); del uvw

        mstable = table(self.msfile, ack=False, readonly=False)
        mstable.putcol('DATA', data)
        mstable.flush()
        mstable.close(); del data

        # ------ Changing Phase Center in FIELD subtable ------ # 
        mstable   = table( os.path.join(self.msfile, 'FIELD'), ack=False, readonly=False)
        delay_tab = mstable.getcol('DELAY_DIR')
        phase_tab = mstable.getcol('PHASE_DIR')
        refer_tab = mstable.getcol('REFERENCE_DIR')
        delay_tab[0, 0, 0] = ra_f # radians !
        delay_tab[0, 0, 1] = dec_f
        phase_tab[0, 0, 0] = ra_f
        phase_tab[0, 0, 1] = dec_f
        refer_tab[0, 0, 0] = ra_f  
        refer_tab[0, 0, 1] = dec_f 
        mstable.putcol('DELAY_DIR', delay_tab)
        mstable.putcol('PHASE_DIR', phase_tab)
        mstable.putcol('REFERENCE_DIR', refer_tab)
        mstable.flush()
        mstable.close()

        # ------ Changing POINTING subtable ------ # 
        mstable   = table( os.path.join(self.msfile, 'POINTING'), ack=False, readonly=False)
        dir_tab   = mstable.getcol('DIRECTION')
        targ_tab  = mstable.getcol('TARGET')
        track_tab = mstable.getcol('TRACKING')
        dir_tab[:,  0, 0] = np.repeat(ra_f,  dir_tab[:,  0, 0].size)
        dir_tab[:,  0, 1] = np.repeat(dec_f, dir_tab[:,  0, 1].size)
        targ_tab[:, 0, 0] = np.repeat(ra_f,  targ_tab[:, 0, 0].size)
        targ_tab[:, 0, 1] = np.repeat(dec_f, targ_tab[:, 0, 1].size)
        track_tab[:] = np.ones(track_tab.size, dtype=bool)
        mstable.putcol('DIRECTION', dir_tab)
        mstable.putcol('TARGET', targ_tab)
        mstable.putcol('TRACKING', track_tab)
        mstable.flush()
        mstable.close()

        self._updateHistory(m='UVW and DATA rephased re-computed to RA='+str(np.degrees(ra_f))+' DEC='+str(np.degrees(dec_f)))
        return

    # --------------------------------------------------------------------------------------------- #
    # ----------------------------------------- Internal ------------------------------------------ #
    def _updateHistory(self, m=''):
        """ Update the history
            m: message
        """
        log = inspect.stack()
        
        mstable = table( os.path.join(self.msfile, 'HISTORY'), ack=False, readonly=True)
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
            mes_table = m
            ori_table = log[1][3]
            pri_table = 'INFO'
            tim_table = Time.now().mjd * 3600. * 24.
        else:
            #cli_table.append( np.array([log[2][4]])  )   # {'shape': [1, 1], 'array': ['        self._makeEmptyMS()\n']}
            app_table.append( log[0][1]  )
            mes_table.append( m )
            ori_table.append( log[1][3] )
            pri_table.append( 'INFO' )
            tim_table = np.append( tim_table, Time.now().mjd * 3600. * 24. )

        mstable = table( os.path.join(self.msfile, 'HISTORY'), ack=False, readonly=False)
        mstable.addrows( nrows=1 )
        mstable.flush()
        mstable.close()

        mstable = table( os.path.join(self.msfile, 'HISTORY'), ack=False, readonly=False)
        #mstable.putcol('CLI_COMMAND', cli_table)
        mstable.putcol('APPLICATION', app_table)
        mstable.putcol('MESSAGE'    , mes_table)
        mstable.putcol('ORIGIN'     , ori_table)
        mstable.putcol('PRIORITY'   , pri_table)
        mstable.putcol('TIME'       , tim_table)
        mstable.flush()
        mstable.close()

        self._taskComplete(task=log[1][3])
        return
    def _azel2radec(self, az, alt, t):
        """ Returns the RA DEC in radians
        """
        tframe  = coord.AltAz(obstime=Time(t, format='jd'), location=self.nancay)
        sAltAz  = coord.SkyCoord(az*u.deg, alt*u.deg, frame=tframe)
        sRaDec  = sAltAz.transform_to(coord.FK5(equinox='J2000'))
        return sRaDec.ra.rad, sRaDec.dec.rad
    def _radecTime(self):
        """ Returns a RA and a DEC arrays, same dimension as the time
            As we expect a fixed 'numeric' pointing / phase center in ALT AZ coords,
            we have to re-compute the EQ coords for the whole obs duration
        """
        # time_list = Time(self.xsttime, format='jd')
        # az_list   = np.repeat(0. , time_list.size)  # np.squeeze(self._azlistana)
        # alt_list  = np.repeat(90., time_list.size) # np.squeeze(self._ellistana)
        # ra, dec   = self._azel2radec(az_list, alt_list, time_list)

        az_list   = np.repeat(0. , self.xsttime.size)  # np.squeeze(self._azlistana)
        alt_list  = np.repeat(90., self.xsttime.size) # np.squeeze(self._ellistana)
        ra, dec   = self._azel2radec(az_list, alt_list, self.xsttime)
        return ra, dec
    def _addCorrected(self):
        """ Copy DATA to CORRECTED_DATA
        """
        if not self.cp_datacol:
            return
        mstable = table(self.msfile, readonly=False)
        data = mstable.getcol('DATA')
        addImagingColumns(self.msfile)
        mstable.putcol('CORRECTED_DATA', data)
        mstable.flush()
        mstable.close()
        return
    def _cleanDirectory(self):
        """ Remove unused sub-items created during the MS construction process
        """
        ant_table = os.path.join(self.savepath, 'ANTENNA')
        try:
            shutil.rmtree(ant_table, ignore_errors=True)
            os.rename( os.path.join(self.savepath, 'makems.cfg'), os.path.join(self.savepath, os.path.basename(self.msfile)[:-3] + '_makems.cfg'))
            os.remove( self.msfile + '.gds' )
            os.remove( self.msfile + '.vds' )
        except:
            print("\t=== WARNING: error trying to clean {} from unused files. ===".format(self.savepath))
        return
    def _taskComplete(self, task):
        """ Display a message to tell everything's allright
        """
        print("\n\t=== Task '{}' completed ({}). ===\n".format(task, Time.now().iso))
        return
    def _splitMS(self, rmms=False):
        """ Split the MS by sub-bands, as LOFAR
            rmms: bool (default = False), remove the Measurement Set once it has correctly been split
        """
        freqt = table( os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=True)
        names = freqt.getcol('NAME')
        freqs = freqt.getcol('REF_FREQUENCY')
        freqt.close()
        desct = table(os.path.join(self.msfile, 'DATA_DESCRIPTION'), ack=False, readonly=True)
        spws  = desct.getcol('SPECTRAL_WINDOW_ID')
        desct.close()

        mstable = table(self.msfile, ack=False, readonly=True)
        for spw, name, freq in zip(spws, names, freqs):
            #output = self.msfile.split('.ms')[0] + '_SB' + str(name) + '_' + str(freq*1.e-6) + 'MHz.ms'
            output = self.msfile.split('.ms')[0] + '_SB' + str(name) + '.ms'
            newms  = mstable.query('DATA_DESC_ID == {}'.format(spw))
            newms.copy(output, deep=True)
            newms.close()
            # Update SPECTRAL_WINDOW and DATA_DESCRIPTION tables to match the frequency selection
            ftable = table( os.path.join(output, 'SPECTRAL_WINDOW'), ack=False, readonly=False)
            ftable.removerows(rownrs=np.arange(self.subbands.size)[np.arange(self.subbands.size) != spw])
            ftable.flush()
            ftable.close()
            mstab  = table(output, ack=False, readonly=False)
            datad  = mstab.getcol('DATA_DESC_ID')
            datad  = np.zeros( datad.shape ) # every element is labelled 'SPW=0' now...
            mstab.putcol('DATA_DESC_ID', datad)
            mstab.flush()
            mstab.close()
            dtable = table( os.path.join(output, 'DATA_DESCRIPTION'), ack=False, readonly=False)
            dtable.removerows(rownrs=np.arange(self.subbands.size)[np.arange(self.subbands.size) != spw])
            spwid = dtable.getcol('SPECTRAL_WINDOW_ID')
            spwid[0] = 0 # set the SPW index to 0 since this is the only frequency
            dtable.putcol('SPECTRAL_WINDOW_ID', spwid)
            dtable.flush()
            dtable.close()
            print("\t=== Measurement set {} written. ===".format(output))
        mstable.close()

        # ------ Remove the MS ------ #
        if rmms:
            sbms = glob.glob( os.path.join(self.savepath, '*_SB*.ms') )
            if len(sbms) != self.subbands.size:
                print("\t=== WARNING: could not find {} sub-band MSs in {}. ===".format(self.subbands.size, self.savepath))
            else:
                try:
                    shutil.rmtree(self.msfile, ignore_errors=True)
                except:
                    print("\t=== WARNING: impossible to remove {} ===".format(self.msfile))
            pass
        return
  
    # --------------------------------------------------------------------------------------------- #
    # ------------------------------------------ Checks ------------------------------------------- #
    def _checkAnt(self):
        """ Check that all required MAs are in the location file
        """
        rfile    = open(self._antlocation, 'r')
        mrinfile = []
        for line in rfile:
            if "DBNAME='NENUFAR'" in line:
                mrinfile.append( int(line.split("'")[1][2:]) )
        for ma in self.miniarrays:
            if ma not in mrinfile:
                raise ValueError("\t=== MA {} not in {} ===".format(ma, self._antlocation))
        return
    def _checkSum(self):
        """ Check that the sum of all auto-correlations imag( XY ) + imag( YX ) is 0
        """
        mstable = table(self.msfile, ack=False, readonly=True)
        mstable = mstable.query('ANTENNA1 == ANTENNA2') # get only auto-correlations
        data = mstable.getcol('DATA')
        mstable.close()
        auto_sum = np.sum(data[:, : , 1].imag + data[:, :, 2 ].imag)
        if auto_sum != 0.:
            raise ValueError("\t=== WARNING: Something's wrong... Failed to verify _checkSum() test! ===")
        else:
            # Everything's alright ;-)
            pass
        return


class NenuMS(object):
    def __init__(self, msfile):
        self.msfile = msfile

    # --------------------------------------------------------------------------------------------- #
    # ------------------------------------------ Plots -------------------------------------------- #
    def plotUV(self, select=None, save=None, plotlambda=False):
        """ Plot the UV plane of the MS
            select: 'DATA_DESC_ID==0' --> SPW 0
        """
        import matplotlib.pyplot as plt

        # ------ Read UV ------ #
        mstable = table(self.msfile, ack=False, readonly=True)
        if select:
            mstable = mstable.query(select)
        uvw = mstable.getcol('UVW')
        if plotlambda:
            # ------ Convert from m to wavelength ------ #
            spwdataid  = mstable.getcol('DATA_DESC_ID')  # Nvis*Ntimes*Kspw
            spwtable   = table(os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=True)
            chanfreq   = spwtable.getcol('CHAN_FREQ')   # Nvis*Ntimes*Kspw x Nchan x 4   in Hz
            desctable  = table( os.path.join(self.msfile, 'DATA_DESCRIPTION'), ack=False, readonly=True)
            spwid      = desctable.getcol('SPECTRAL_WINDOW_ID')
            frequency  = np.take(chanfreq, spwid[spwdataid], axis=0)
            wavelength = const.c.value / frequency
            uvw *= 1./wavelength
        else:
            pass
        u   = np.append( -uvw[:, 0], uvw[:, 0])
        v   = np.append( -uvw[:, 1], uvw[:, 1])
        mask = (u != 0.) & (v != 0) # no auto-correlations
        xmin, xmax = u.min(), u.max()
        ymin, ymax = v.min(), v.max() 
        
        # ------ Make a density plot ------ #
        fig, ax = plt.subplots()
        hb = ax.hexbin(u[mask], v[mask], gridsize=400, mincnt=1, bins='log', cmap='Blues', edgecolor='none')
        ax.margins(0) # must have!
        ax.axis([xmin, xmax, ymin, ymax])
        ax.set_title('UV Coverage')
        if plotlambda:
            ax.set_xlabel('u ($\\lambda$)')
            ax.set_ylabel('v ($\\lambda$)')
        else:
            ax.set_xlabel('u (m)')
            ax.set_ylabel('v (m)')
        cb = fig.colorbar(hb, ax=ax)
        cb.set_label('log$_{10}$(density)')
        if save:
            plt.savefig(save)
        else:
            plt.show()
        return
    def plotLines(self, select=None, colname='CORRECTED_DATA'):
        """ Cover of drawMS by Cyril Tasse - adapted for NenuFAR MS
            import nenupyms; a=nenupyms.NenuMS('NENUFAR_XST_20180603_190400.ms'); a.plotLines(select='DATA_DESC_ID==0')
        """
        self.colname = colname
        self.select  = select  
        self._readMS()
        if self.n_channel > 1:
            print("\t=== Multiple frequencies detected in {}, please use select='DATA_DESC_ID==0' for e.g. ===".format(self.msfile))
            return

        # ------ Find Sources ------ #
        self.ntheta = 20
        self.dt     = 500 # seconds
        self.npeaks = 7 # dunno...
        self.fillf  = 1. # fill_factor
        self.time_range_min = 0.
        self.time_range_max = 1.
        self.wmax = 1000000.
        self.blcut_min = 10#100
        self.blcut_max = 10000000#100000
        self.mask_freq = 5
        self.snrcut = 5.
        self.freq_cut = -1
        
        # Catalog initialization
        dtype_line=[('x', np.float32, (self.ntheta,)),
            ('y',np.float32,(self.ntheta,)),
            ('r',np.float32,(self.ntheta,)),
            ('theta',np.float32,(self.ntheta,))]
        self.cat_lines = np.zeros( (int(3.*self.n_baseline*(self.indt.size/self.dt +1.)*self.npeaks),),
            dtype=[('bln', np.int32), ('timeslot',np.float), ('line', dtype_line), ('freq', np.float32), ('flag', np.bool)])
        self.cat_lines = self.cat_lines.view(np.recarray)
        self.cat_lines.timeslot = -1.
        self.cat_lines.line.x   = -1.
        self.cat_lines.line.y   = -1.
        self.cat_lines.flag     = False

        # Calculation
        ss = range(0, self.indt.size, self.dt)
        if ss[-1] != (self.indt.size-1):
            ss.append(self.indt.size-1)
        if len(ss) < 2:
            print("\t==== Time bin too long, there was no selected time windows ===")
            return

        coutnbl   = 0
        index_cat = 0
        bar = ChargingBar('Calculating lines....', max=self.n_baseline)
        for bl in xrange(self.n_baseline):
            for i in xrange( len(ss)-1 ):
                if np.random.rand(1)[0] > self.fillf:
                    continue
                x, y, flist, rout, theta_rad, snr = self._testFind(bl, ss[i], ss[i+1])
                if x.shape != (1, 1):
                    coutnbl += 1
                    for ii in xrange(len(x)):
                        for iii in xrange(2):
                            line = 0
                            cont = True
                            while cont:
                                rr = np.sqrt(x[ii][iii]**2 + y[ii][iii]**2)
                                mask, cont = self._giveMask(np.abs(rr - 0.555) > 0.0001, line)
                                line += 1
                                if cont:
                                    self.cat_lines.line.x[index_cat][mask] = np.float32(x[ii][iii][mask])
                                    self.cat_lines.line.y[index_cat][mask] = np.float32(y[ii][iii][mask])
                                    self.cat_lines.line.theta  [index_cat][mask] = np.float32(theta_rad[mask])
                                    self.cat_lines.line.r      [index_cat][mask] = np.float32(rout[ii][iii][mask])
                                    self.cat_lines.freq        [index_cat]       = flist[ii]
                                    self.cat_lines.bln         [index_cat]       = bl
                                    self.cat_lines.timeslot    [index_cat]       = i
                                    index_cat += 1
            bar.next()
        bar.finish()
        cat = self.cat_lines[ self.cat_lines.timeslot != -1. ]
        np.save('catalog_lines', cat)
        del (self.cat_lines, cat)

        # ------ Plot ------ #        
        self._plotGridded()
        return

    # --------------------------------------------------------------------------------------------- #
    # ----------------------------------------- Internal ------------------------------------------ #
    def _readMS(self):
        """ Read the MS
        """
        print("\t=== Reading {}...===".format(self.msfile))
        mstable = table(os.path.join(self.msfile, 'ANTENNA'), ack=False, readonly=True)
        if self.select:
            try: mstable = mstable.query(self.select)
            except: pass
        self.n_antenna  = mstable.getcol('POSITION').shape[0]
        self.n_baseline = (self.n_antenna*(self.n_antenna-1))/2 + self.n_antenna
        mstable.close()
        mstable = table( os.path.join(self.msfile, 'SPECTRAL_WINDOW'), ack=False, readonly=True)
        if self.select:
            try: mstable = mstable.query(self.select)
            except: pass
        self.ref_freq   = mstable.getcol('REF_FREQUENCY')[0] # Hz
        self.wavelength = const.c.value / self.ref_freq 
        self.wavelength_chan = const.c.value / mstable.getcol('CHAN_FREQ')[0]
        self.n_channel  = self.wavelength_chan.shape[0]
        mstable.close()
        mstable = table( os.path.join(self.msfile, 'FIELD'), ack=False, readonly=True)
        if self.select:
            try: mstable = mstable.query(self.select)
            except: pass
        self.rarad, self.decrad = mstable.getcol('PHASE_DIR')[0][0]
        mstable.close()
        self.radeg, self.decdeg = np.degrees(self.rarad), np.degrees(self.decrad)
        mstable = table(self.msfile, ack=False, readonly=True)
        if self.select:
            try: mstable = mstable.query(self.select)
            except: pass
        self.vis  = mstable.getcol(self.colname)
        self.uvw  = mstable.getcol('UVW') / self.wavelength_chan[-1]
        self.ant0 = mstable.getcol('ANTENNA1')
        self.ant1 = mstable.getcol('ANTENNA2')
        self.time = mstable.getcol('TIME')
        self.indt = self.time[ 0::self.n_baseline ]
        self.flag = mstable.getcol('FLAG')
        self.vis[self.flag] = 0.
        return
    def _radec2lm(self, ra, dec):
        l = np.cos(dec) * np.sin(ra - self.rarad)
        m = np.sin(dec) * np.cos(self.decrad) - np.cos(dec) * np.sin(self.decrad) * np.cos(ra - self.rarad)
        return l, m
    def _testFind(self, bl, s0, s1):
        """ Need in plotLines()
        """
        import scipy.fftpack
        import scipy.interpolate

        uvw   = self.uvw[  bl::self.n_baseline ]
        time  = self.time[ bl::self.n_baseline ]
        vis   = self.vis[  bl::self.n_baseline, :, 0] + self.vis[ bl::self.n_baseline, :, 3]
        wmean = np.abs(np.mean(uvw[:, 2]))

        # Baselines
        baseline_dist = np.sqrt(uvw[:, 1]**2 + uvw[:, 0]**2)
        baseline_dist_mean = np.mean(baseline_dist)
        ant0 = self.ant0[bl]
        ant1 = self.ant1[bl]
        tbin = (s0 + s1)/2.
        condtime = (tbin < self.time_range_min*uvw.shape[0]) | (tbin > self.time_range_max*uvw.shape[0])
       
        if (wmean>self.wmax)|(uvw[0, 0]==0.)|(baseline_dist_mean<self.blcut_min)|(baseline_dist_mean>self.blcut_max)|condtime:
            return np.array([[0.]]), np.array([[0]]), [0], [0], [0], [0]

        uvw  = uvw[s0:s1]
        time = time[s0:s1]
        vsel = np.array(vis[s0:s1], dtype=np.complex)
        vsel[np.isnan(vsel)] = 0.

        # Fourier
        fft    = scipy.fftpack.fftshift(scipy.fftpack.fft(vsel[:, 0]))
        ff0    = self.wavelength_chan[0]*scipy.fftpack.fftshift(scipy.fftpack.fftfreq(vsel.shape[0], d=time[1]-time[0]))
        # fft_ch = scipy.fftpack.fftshift(scipy.fftpack.fft(vsel[:,0]))
        # ff_ch  = self.wavelength_chan[0]*scipy.fftpack.fftshift(scipy.fftpack.fftfreq(vsel.shape[0],d=time[1]-time[0]))
        # inter  = scipy.interpolate.interp1d(ff_ch, fft_ch)
        # fft += inter(ff0)
        absfft = np.abs(fft)
        maxff  = np.max(absfft)
        stdff  = np.std(absfft)
        meanff = np.median(absfft)
        ff     = ff0 / self.wavelength_chan[0]

        u     = uvw[:,0]
        v     = uvw[:,1]
        w     = uvw[:,2]
        du    = u[0:-1] - u[1:]
        dv    = v[0:-1] - v[1:]
        dw    = w[0:-1] - w[1:]
        dt    = time[0:-1] - time[1:]
        du_dt = np.median(du / dt)
        dv_dt = np.median(dv / dt)
        dw_dt = np.median(dw / dt)
        
        theta_rad = np.sort( np.array(range(self.ntheta))/(self.ntheta-1.)*2.*np.pi )
        
        if (maxff < self.snrcut*stdff)|(maxff==0.)|(stdff < 1.e-6):
            return np.array([[0.]]), np.array([[0]]), [0], [0], [0], [0]

        acoef  = -0.5 * dw_dt
        bcoef  = ( du_dt*np.cos(theta_rad) + dv_dt*np.sin(theta_rad) )
        npeaks = 0
        xout   = []
        yout   = []
        flist  = []
        rout   = []
        snr    = []
        while True:
            theta_rad +=+ np.random.rand(1)[0]*2.*np.pi
            npeaks += 1
            if npeaks > self.npeaks: break
            if maxff < meanff + self.snrcut*stdff: break
            ind  = np.argmax(absfft)
            freq = ff[ind]
            # absfft[ind - self.mask_freq : ind + self.mask_freq] = 0.
            if np.abs(freq) < self.freq_cut: continue
            ccoef = -freq

            rs_sol = np.zeros((2, theta_rad.shape[0]), dtype=np.float) - 0.555
            if acoef != 0.:
                D = bcoef**2 - 4.*acoef*ccoef
                ind0 = np.where(D>0.)[0]
                if ind0.shape[0]>0:
                    rs_sol[0, ind0] = -(bcoef[ind0] - np.sqrt(D[ind0]))/(2.*acoef)
                    rs_sol[1, ind0] = -(bcoef[ind0] + np.sqrt(D[ind0]))/(2.*acoef)
            else:
                rs_sol[0, :]= -ccoef/bcoef

            for sol in xrange(2):
                ind = np.where((np.abs(rs_sol[sol, :]) > -0.3) & (np.abs(rs_sol[sol, :]+0.555) > 0.00001))[0]
                replace = 0
                if len(ind) > 0:
                    rs_sol[sol, ind] = -0.555
                    rs = self._solveRadius(du_dt, dv_dt, dw_dt, freq, theta_rad[ind])
                    ind0 = np.where(rs[0, :] != -0.555)[0]
                    if len(ind0) > 0:
                        rs_sol[0, ind] = rs[0,:]
                        replace=1
                    ind0 = np.where(rs[1, :] != -0.555)[0]
                    if len(ind0) > 0:
                        rs_sol[1, ind]=rs[1,:]
                        replace=2
                
            ind = np.where(rs_sol > 1.)[0]
            if len(ind) > 0: stop
 
            xout.append(  rs_sol*np.cos(theta_rad) )
            yout.append(  rs_sol*np.sin(theta_rad) )
            rout.append(  rs_sol )
            flist.append( freq )
            snr.append(   maxff/stdff )                        
        if xout == 0: 
            return np.array([[0.]]), np.array([[0]]), [0], [0], [0], [0]

        return np.array(xout), np.array(yout), flist, rout, theta_rad, snr
    def _giveMask(self, mask_bool, n):
                indn = -1
                ind  = 0
                aold = False
                mask_out=np.zeros(mask_bool.shape,dtype=np.bool)
                while (mask_bool[ind] != True) | (indn != n):
                    if (mask_bool[ind] == True) & (aold == False):
                        indn += 1
                    if ind == mask_bool.shape[0]-1:
                        return mask_out, False
                    aold = mask_bool[ind]
                    ind += 1
                ind -= 1
                while mask_bool[ind]:
                    mask_out[ind] = True
                    if ind==mask_bool.shape[0]-1: return mask_out, True
                    ind += 1
                return mask_out, True
    def _solveRadius(self, du_dt, dv_dt, dw_dt, freq, theta):
        import scipy.optimize
        rads = np.linspace(0., 0.999, 10)
        rout = np.zeros((2, theta.shape[0]), dtype=np.float) - 0.555
        for i in range(theta.shape[0]):
            def freq_lmn(rad, du_dt, dv_dt, dw_dt, freq, theta): 
                return rad*(du_dt*np.cos(theta)+dv_dt*np.sin(theta))+dw_dt*(np.sqrt(1.-rad**2)-1.)-freq
            yy   = rads*(du_dt*np.cos(theta[i]) + dv_dt*np.sin(theta[i])) + dw_dt*(np.sqrt(1. - rads**2) - 1.) - freq
            yyp  = du_dt*np.cos(theta[i]) + dv_dt*np.sin(theta[i]) - dw_dt*rads/np.sqrt(1. - rads**2) # derivative
            cond = [yyp < 0., yyp > 0.]
            for ii in xrange(2):
                ind, end = self._giveMask(cond[ii],0)
                if not end: continue
                radssel = rads[ind]
                yysel   = yy[ind]
                if (np.min(yysel) < 0.) & (np.max(yysel) > 0.):
                    rout[ii, i] = scipy.optimize.brentq(freq_lmn, np.min(radssel), np.max(radssel), 
                        args=(du_dt, dv_dt, dw_dt, freq, theta[i]), xtol=1e-3, rtol=1e-3, maxiter=100, full_output=False, disp=True)
        ind = np.where(rout == -0.555)
        if freq==0.:
            rout[ind]=0.
        return rout
    def _plotGridded(self):
        """ Plot the lines
        """
        import matplotlib.pyplot as plt
        cat = np.load('catalog_lines.npy')
        cat = cat.view(np.recarray)

        npix = 2000
        img  = np.ones((npix,npix))
        n_interp_line = 100
        dens_interp_line = 1e-4
        consize = 2

        falsecount = 0
        for c in range(cat.shape[0]):
            ind = np.where(
                (cat.line.x[c] != -1.) &
                (cat.line.x[c] != np.nan) &
                (cat.line.y[c] != np.nan) &
                (cat.line.x[c] != 1.) &
                (cat.line.y[c] != -1.) &
                (cat.line.y[c] != 1.))
            if ind[0].shape[0] > 0:
                for i in xrange(ind[0].shape[0] - 1):
                    xm = cat.line.x[c][ind[0]][i:i+2]
                    ym = cat.line.y[c][ind[0]][i:i+2]
                    dist = np.sqrt((xm[0]-xm[1])**2 + (ym[0]-ym[1])**2)
                    indu = xm.argsort()
                    xm = xm[indu]
                    ym = ym[indu]
                    if dist==0: 
                        falsecount += 1
                        continue
                    n_interp_line = dist/dens_interp_line
                    xs = np.arange(np.min(xm), np.max(xm), (np.max(xm)-np.min(xm))/(n_interp_line))
                    ys = np.interp(xs, xm, ym)
                    x = npix*(xs+1.)/2.
                    y = npix*(ys+1.)/2.
                    indnan = np.where((x!=np.nan) & (y!=np.nan))[0]
                    img[np.int16(x[indnan]), np.int16(y[indnan])] += 1
        img = np.log10(img)
        plt.imshow(img.T, origin='lower', extent=(-1., 1., -1., 1.), cmap='bone_r', interpolation='gaussian')
        plt.xlabel("L")
        plt.ylabel("M")
        plt.title("Image plane")
        theta = np.arange(0., 2.*np.pi, 0.01)
        plt.plot(np.cos(theta), np.sin(theta), color='black')

        # Overplot A team positions
        ateam = [ {'name' : 'Cas A', 'ra' : 6.123487680622104,  'dec' : 1.0265153995604648, 'tag': 0},
              {'name' : 'Cyg A', 'ra' : 5.233686575770755,  'dec' : 0.7109409582180791, 'tag': 0},
              {'name' : 'Tau A', 'ra' : 1.4596748493730913, 'dec' : 0.38422502335921294, 'tag': 0},
              {'name' : 'Her A', 'ra' : 4.4119087330382163, 'dec' : 0.087135562905816893, 'tag': 0},
              {'name' : 'Vir A', 'ra' : 3.276086511413598,  'dec' : 0.21626589533567378, 'tag': 0}]
        import pyrap.quanta as qa
        import pyrap.measures as pm
        import ephem
        me = pm.measures()
        dict_time_start_mjd = me.epoch('utc', qa.quantity(self.time[0], 's'))
        time_start_mjd=dict_time_start_mjd['m0']['value']
        jd=time_start_mjd + 2400000.5 - 2415020
        d=ephem.Date(jd)
        sun, jup = ephem.Sun(), ephem.Jupiter()
        sun.compute(d)
        jup.compute(d)
        ateam.append({'name': 'Sun', 'ra': sun.ra, 'dec': sun.dec, 'tag': 1})
        ateam.append({'name': 'Jupiter', 'ra': jup.ra, 'dec': jup.dec, 'tag': 1})

        for i in xrange(len(ateam)):
            pi2 = np.pi/2.
            ra  = ateam[i]['ra']
            dec = ateam[i]['dec']
            cosd = np.cos(pi2-dec)*np.cos(pi2-self.decrad) + np.sin(pi2-dec)*np.sin(pi2-self.decrad)*np.cos(ra-self.rarad)
            d = np.arccos(cosd)
            if d < np.pi/2.:
                l, m = self._radec2lm(ra, dec)
                name = ateam[i]['name']
                plt.plot([-l], [m], marker='o', color='#d62728', fillstyle='none')
                #plt.scatter([-l], [m], facecolors='none', edgecolors='#d62728')
                plt.text(-l, m, ateam[i]['name'], color='#d62728')

        plt.savefig('Lines.pdf')
        return





