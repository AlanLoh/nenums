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
    raise ImportError("\n\t=== Unable to import the pyrap module ===")


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
        return

