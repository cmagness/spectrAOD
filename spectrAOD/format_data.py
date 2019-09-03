#!/usr/bin/env python

"""This module is for formatting input data into an object that stores
primarily the wavelength and flux data, as well
as any relevant target information. This object is then used by
measure_aod.py to perform the measurements."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import os
import sys
import glob

import numpy as np
from astropy import constants
from astropy.io import fits
from astropy.io import ascii

from . import SETTINGS
from .spectrum_classes import X1DSpectrum, ASCIISpectrum

INPUTS = SETTINGS["inputs"]
DATADIR = INPUTS["datadir"]
PARAMETERS = SETTINGS["parameters"]

C = constants.c.to('km/s').value  # km/s

N = 3.768e14  # proportionality constant -> (m_e * c)/(pi * e**2)

# at some point need to check and see if outdir (and datadir, really too)
# exist. if not, create outdir
# add logger instead of print statements

# notes on class functionality:
# add functionality to build fits files as a method of class objects


# --------------------------------------------------------------------------- #


def build_spectrum(datadir=DATADIR, ins=PARAMETERS["instrument"],
                   file=PARAMETERS["filetype"], grating=PARAMETERS["grating"]):
    # this function should build and return the appropriate Spectrum object
    # for measure_aod.py
    # default inputs for instrument and file type are COS and X1D at the
    # moment, for testing

    spectrum = None

    # this might (?) need handling for other inputs other than just COS as
    # well. at some point.
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
                        wave = np.array(target_data[
                                            "WAVELENGTH"].ravel())
                        # best way to do this, ravelling here?
                        flux = np.array(target_data["FLUX"].ravel())
                        error = np.array(target_data["ERROR"].ravel())
                        spectrum = X1DSpectrum(target, wave, flux, error)
        elif file == "BART":
            search_string = "_spec-{}".format(grating)
            asciis = glob.glob(datadir + "*" + search_string)
            # multiple asciis??
            for ascii_file in asciis:
                basename = os.path.basename(ascii_file)
                # this can be done better
                target = basename.replace(search_string, "")
                data = ascii.read(ascii_file, names=["wave", "flux", "error"])
                if len(data.columns) > 3:
                    raise IndexError("Too many columns")
                    # choose a better error type?
                wave = np.array(data.columns["wave"])
                flux = np.array(data.columns["flux"])
                error = np.array(data.columns["error"])
                spectrum = ASCIISpectrum(target, wave, flux, error)
        elif file == "BART-N":
            search_string = "_spec-{}-N".format(grating)
            asciis = glob.glob(datadir + "*" + search_string)
            for ascii_file in asciis:
                basename = os.path.basename(ascii_file)
                # this can be done better
                target = basename.replace(search_string, "")
                data = ascii.read(ascii_file, names=["wave", "flux", "error"])
                if len(data.columns) > 3:
                    raise IndexError("Too many columns")
                    # choose a better error type?
                wave = np.array(data.columns["wave"])
                flux = np.array(data.columns["flux"])
                error = np.array(data.columns["error"])
                spectrum = ASCIISpectrum(target, wave, flux, error)
        else:
            print("Other file types are not yet supported at this time.")
    else:
        print("This Instrument is not yet supported at this time.")

    # this needs to be addressed with proper error handling so that it is not
    # referenced before assignment
    if not spectrum:
        sys.exit("Spectrum object not built. Unable to proceed. Exiting...")

    return spectrum


# --------------------------------------------------------------------------- #
