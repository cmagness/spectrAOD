#!/usr/bin/env python

"""This module is for formatting input data into an object that stores primarily the wavelength and flux data, as well
as any relevant target information. This object is then used by measure_aod.py to perform the measurements."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"


import glob
from astropy.io import fits

DATADIR = "/user/cmagness/fermi/data/15339/UVQSJ191928-295808/x1d/"
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

    def __init__(self, target, wave, flux):
        self.target = target
        self.wave = wave
        self.flux = flux

    def tofits(self):
        # this method should generate a fits file
        pass

# ----------------------------------------------------------------------------------------------------------------------


class X1DSpectrum(BaseSpectrum):

    def __init__(self, *args):
        super().__init__(*args)

    # we're gonna wanna add other things we can capture with x1d spectra here later. which will also require updates to
    # that part of the if statement in format data


# ----------------------------------------------------------------------------------------------------------------------

# add classes for other data types

# ----------------------------------------------------------------------------------------------------------------------


def format_data(datadir=DATADIR, ins="COS", file="X1DSUM", grating="G130M"):
    # this function should build and return the appropriate Spectrum object for measure_aod.py
    # default inputs for instrument and file type are COS and X1D at the moment, for testing

    # this might (?) need handling for other inputs other than just COS as well. at some point.
    if ins == "COS":
        if file == "X1DSUM":
            x1dsums = glob.glob(datadir + "*x1dsum.fits")
            for x1dsum in x1dsums:
                with fits.open(x1dsum) as f:
                    prhd = f["PRIMARY"].header
                    opt_elem = prhd["OPT_ELEM"]
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