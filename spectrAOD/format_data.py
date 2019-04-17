#!/usr/bin/env python

"""This module is for formatting input data into an object that stores primarily the wavelength and flux data, as well
as any relevant target information. This object is then used by measure_aod.py to perform the measurements."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import glob
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

DATADIR = "/user/cmagness/fermi/data/15339/UVQSJ191928-295808/x1d/"

C = 2.99792458e5  # km/s


# OUTDIR = "/user/cmagness/fermi/data/out/" this is unused at this time but will be eventually

# at some point need to check and see if outdir (and datadir, really too) exist. if not, create outdir
# add logger instead of print statements

# notes on class functionality:
# should have attributes for wavelength, flux at a minimum
# also would be great to have target name, l, b coordinates. or conversion from ra & dec if that information is
# available
# add functionality to build fits files as a method of class objects


# ----------------------------------------------------------------------------------------------------------------------


class BaseSpectrum:
    """This is the base class for the Spectrum objects that get used to store data for this package."""

    # NEED TO ADD SOME ERROR HANDLING IF METHODS ARE CALLED BEFORE ATTRIBUTES USED ARE DEFINED

    def __init__(self, target, wave, flux):
        self.target = target
        self.wave = wave
        self.flux = flux
        self.doublet = False
        self.velocity = None
        self.ra = None
        self.dec = None
        self._skycoords = None
        # ADD A FLAG FOR INSIDE VS OUTSIDE FERMI BUBBLES

    @property
    def l(self):
        return self._skycoords.galactic.l.value

    @property
    def b(self):
        return self._skycoords.galactic.b.value

    def set_doublet(self, doublet):
        self.doublet = doublet

    def calculate_velocity(self, ion_wv: dict):
        # this method calculates the wavelength in velocity space wrt to the ion of interest and stores it
        for wv in ion_wv:
            z_array = (self.wave - wv) / wv
            vel_array = C * z_array
            self.velocity.append(vel_array)
        # eventually will need to keep track of ion associated with the velocity array
        # for now it is fine because they are associated with the same ion in the doublet case

    def get_coords(self, target_list):
        # this method should get the RA & DEC from the coordinates list
        df_targets = pd.read_csv(target_list, index_col=0)
        # need to add some error handling for if target name is not found in list
        mask = df_targets["Target"] == self.target
        self.ra = df_targets.loc[mask]["RA"]
        self.dec = df_targets.loc[mask]["DEC"]
        self._skycoords = SkyCoord(ra=self.ra * u.degree, dec=self.dec * u.degree)

    def lsr_correct_velocity(self):
        # this method should apply the lsr correction to the heliocentric coordinates
        l_radians = self.l * np.pi / 180.0
        b_radians = self.b * np.pi / 180.0
        vel_corr = (9.0 * np.cos(l_radians) * np.cos(b_radians)) + (12.0 * np.sin(l_radians) * np.cos(b_radians)) + \
                   (7.0 * np.sin(b_radians))
        self.velocity = self.velocity + vel_corr

    def to_fits(self):
        # this method should generate a fits file
        pass


# ----------------------------------------------------------------------------------------------------------------------


class X1DSpectrum(BaseSpectrum):
    """This class inherits the base Spectrum class for specifically x1dsum files and will have methods for holding
    other x1d specific information."""

    def __init__(self, *args):
        super().__init__(*args)

    def x1d_specs(self):
        # this method could hold other x1d specific stuff, like header information of interest
        pass

    # we're gonna wanna add other things we can capture with x1d spectra here later. which will also require updates to
    # that part of the if statement in format data


# ----------------------------------------------------------------------------------------------------------------------

# add classes for other data types

# ----------------------------------------------------------------------------------------------------------------------


def format_data(datadir=DATADIR, ins="COS", file="X1DSUM", grating="G130M"):
    # this function should build and return the appropriate Spectrum object for measure_aod.py
    # default inputs for instrument and file type are COS and X1D at the moment, for testing

    spectrum = []

    # this might (?) need handling for other inputs other than just COS as well. at some point.
    if ins == "COS":
        if file == "X1DSUM":
            x1dsums = glob.glob(datadir + "*x1dsum.fits")
            for x1dsum in x1dsums:
                with fits.open(x1dsum) as f:
                    prhd = f["PRIMARY"].header
                    opt_elem = prhd["OPT_ELEM"]
                    # WHAT IF IT FINDS MULTIPLE X1DSUMS HERE?
                    if opt_elem == grating:
                        target_data = f["SCI"].data
                        target = prhd["TARGNAME"]
                        wave = target_data["WAVELENGTH"]
                        flux = target_data["FLUX"]
                        spectrum = X1DSpectrum(target, wave, flux)
        else:
            print("Other file types are not yet supported at this time.")
    else:
        print("This Instrument is not yet supported at this time.")

    # this needs to be addressed with proper error handling so that it is not referenced before assignment
    if not spectrum:
        spectrum = []
        print("Spectrum object not built.")

    return spectrum
