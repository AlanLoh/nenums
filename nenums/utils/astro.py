#! /usr/bin/python3.6
# -*- coding: utf-8 -*-

"""
Functions to do astronomical stuff
        by A. Loh
"""

from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy import coordinates as coord
from astropy import constants as const

__author__ = ['Alan Loh']
__copyright__ = 'Copyright 2018, nenums'
__credits__ = ['Alan Loh']
__license__ = 'MIT'
__version__ = '0.0.1'
__maintainer__ = 'Alan Loh'
__email__ = 'alan.loh@obspm.fr'
__status__ = 'WIP'
__all__ = ['nancay', 'altaz2radec']


def nancay():
    """ NenuFAR's position
        
        Returns
        -------
        Coordinate object 
    """
    return coord.EarthLocation(lat=47.376511*u.deg, lon=2.1924002*u.deg)


def lightSpeed():
    """ Get the light speed

        Returns
        -------
        * **c** : float
            Speed of light in m/s
    """
    return const.c.value


def altaz2radec(az, alt, t):
    """ Convert Azimuth-Elevation coordinates into RA-Dec equatorial coordinates
        It is expected that the observer is NenuFAR in Nancay
        The conversion is computed for a sky set at time t
        
        Parameters
        ----------
        az : float
            Azimuth in degrees
        alt : float
            Elevation in degrees
        t : str or Time
            Time at which to compute the corresponding RA Dec

        Returns
        -------
        RA : float
            Right Ascension in radians
        Dec : float
            Declination in radians
    """
    if not isinstance(t, Time):
        try:
            t = Time(t)
        except:
            raise ValueError("\n\t=== Time syntax should be 'YYYY-MM-DD hh:mm:ss' ===") 
    frame = coord.AltAz(obstime=t, location=nancay())
    altaz = coord.SkyCoord(az*u.deg, alt*u.deg, frame=frame)
    radec = altaz.transform_to(coord.FK5(equinox='J2000'))
    return radec.ra.rad, radec.dec.rad






