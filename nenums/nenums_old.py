#!/usr/bin/env python
__author__ = "Alan Loh"
__credits__ = ["Alan Loh"]
__version__ = "1.0.1"
__maintainer__ = "Alan Loh"
__email__ = "alan.loh@obspm.fr"
__status__ = "Production"

"""
=========================================================================
                                DESCRIPTION
    Transform a XST file into a Measurement Set
    It also recomputes UVW coorfinates and numerically re-phase the visibility at a new phase center (ra, dec)

    Arguments:
        - obs   (needed)   : original XST file (string)
        - store (needed)   : path (existing or not) where the new Measurement Set will be created
        - track (optionnal): recompute or not the coordinates and rephase the visibilities
        - ra    (optionnal): new phase center RIGHT ASCENSION (in degrees)
        - dec   (optionnal): new phase center DECLINATION (in degrees)
        - split (optionnal): split the ms by sub-bands 
    
    Example:
    nenums.py --obs /usr/data/xst_observation.fits --store /usr/path/to/save --ra 300. --dec 45.
    Or (short version):
    nenums.py -o /usr/data/xst_observation.fits -s /usr/path/to/save -r 300. -d 45.
=========================================================================
"""

import nenupyms
import argparse

# =========================================================================
# =========================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--obs',   type=str,   required=True,                help='(REQUIRED) original XST file (string)')
    parser.add_argument('-s', '--store', type=str,   required=True,                help='(REQUIRED) path (existing or not) where the new Measurement Set will be created')
    parser.add_argument('-t', '--track', type=int,  required=False, default=1, help='(OPTION, default=1) recompute or not the coordinates and rephase the visibilities')
    parser.add_argument('-r', '--ra',    type=float, required=False, default=None, help='(OPTION, default=None --> RA read in XST.fits) new phase center RIGHT ASCENSION (in degrees)')
    parser.add_argument('-d', '--dec',   type=float, required=False, default=None, help='(OPTION, default=None --> DEC read in XST.fits) new phase center DECLINATION (in degrees)')
    parser.add_argument('-sp', '--split',   type=int, required=False, default=1, help='(OPTION, default=1) Split the MS by Sub-Bands as LOFAR')
    parser.add_argument('-c', '--corrected',   type=int, required=False, default=0, help='(OPTION, default=0) Add a CORRECTED_DATA column (identical to DATA)')
    args = parser.parse_args()

    obs = nenupyms.NenuFits(obsfile=args.obs, savepath=args.store)
    if (args.ra is not None) and (args.dec is not None):
        obs.forceradec = (args.ra, args.dec)
    #if not (args.track):
    #    obs.do_tracking = False
    #if not (args.corrected):
    #    obs.cp_datacol = False
    obs.do_tracking = bool(args.track)
    obs.cp_datacol = bool(args.corrected)
    obs.createMs(sbsplit=bool(args.split))
    del obs