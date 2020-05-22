#!/usr/bin/env python

"""This module is for measuring the apparent optical depth of absorption lines
in COS Fermi Bubble spectra. The method of measurement is outlined in Savage &
Sembach 1991. Eventually, this will be extended to measure apparent optical
depth of other spectra as well."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import sys
import argparse
import logging

from . import SETTINGS
from .format_data import build_spectrum
from .spectrum_classes import Helper

INPUTS = SETTINGS["inputs"]
DATADIR = INPUTS["datadir"]
TARGETS = INPUTS["targets"]
PARAMETERS = SETTINGS["parameters"]
DEFAULTS = SETTINGS["defaults"]

# set up log file. will overwrite with each new run
# log file has DEBUG level + written to it, package, module, function
# console has only INFO level + written, basic format
logging.basicConfig(filename="logfile.log",
                    format="%(levelname)-4s[%(module)s.%(funcName)s.%("
                           "lineno)d]:%(message)s", filemode="w",
                    level=logging.DEBUG)
LOGGER = logging.getLogger(__name__)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
LOGGER.addHandler(console)


# --------------------------------------------------------------------------- #


def main():
    # argparse
    args, spectrum = parse()
    LOGGER.info("Spectrum object successfully built.")
    # LSR correction
    spectrum = lsr_correct(args, spectrum)
    LOGGER.info("Spectrum LSR corrected.")
    # continuum fit
    spectrum = continuum_fit(spectrum)
    LOGGER.info("Continuum fit calculated.")
    # measure aod/acd/ew
    # set measurements back in spectrum object from helper object
    spectrum = measure(args, spectrum)
    # generate table
    spectrum.generate_table(args.vel_min, args.vel_max)
    LOGGER.info("spectrAOD complete.")
    return 0


# --------------------------------------------------------------------------- #


def parse():
    """This function is for parsing the inputs of interest."""

    parser = argparse.ArgumentParser(description="spectrAOD")

    # for debug only
    # parser.add_argument("--ion", default="SiIV", type=str,
    #                     help="absorption line feature of interest")

    parser.add_argument("ion", type=str,
                        help="absorption line feature of interest")

    parser.add_argument("--instrument", default=PARAMETERS["instrument"],
                        type=str, help="observational instrument data is taken"
                                       "on")
    parser.add_argument("--filetype", default=PARAMETERS["filetype"], type=str,
                        help="filetype of data, takes X1DSUM (for x1ds or "
                             "x1dsums), BART, and BART-N")
    parser.add_argument("--vel_min", default=PARAMETERS["vel_min"], type=int,
                        help="velocity minimum in window of interest around "
                             "ion, in km/s")
    parser.add_argument("--vel_max", default=PARAMETERS["vel_max"], type=int,
                        help="velocity maximum in window of interest around "
                             "ion, in km/s")
    parser.add_argument("--grating", default=PARAMETERS["grating"], type=str,
                        help="grating of observation mode")
    parser.add_argument("--redshift", default=PARAMETERS["redshift"],
                        type=float, help="redshift of wavelength "
                                         "region/target")
    # need to do something about the grating bc not every dataset will have
    # a grating in the future

    # need to set up better debugging system than defaults. debugging flag?
    # need to set default grating to something for creating Spectrum object
    # NEED TO DO SOMETHING ABOUT THE FACT THAT COULD HAVE ION LIST IN THE
    # FUTURE

    # grating evaluation needs some work because it will have to be checked
    # that the ion of interest falls in wavelength
    # range of the grating

    args = parser.parse_args()
    LOGGER.info("Initialized Arguments: \n"
                "Ion: {} \n"
                "Instrument: {} \n"
                "Filetype: {} \n"
                "Velocity Window: [{}, {}] \n"
                "Grating: {} \n"
                "Redshift: {}".format(args.ion, args.instrument.upper(),
                                      args.filetype.upper(), args.vel_min,
                                      args.vel_max, args.grating.upper(),
                                      args.redshift))
    LOGGER.warning("If these are *not* the parameters you intended to use, "
                   "please perform a new measurement. Future versions of "
                   "this package will support mid process adjustments of "
                   "these parameters.")
    spectrum = build_spectrum(DATADIR, args.instrument.upper(),
                              args.filetype.upper(), args.grating.upper(),
                              args.redshift)

    return args, spectrum


# --------------------------------------------------------------------------- #


def lsr_correct(args, spectrum):
    """This function performs the lsr correction"""

    # find ion wavelength from ions.csv
    spectrum.get_ions(args.ion)
    # transform wavelength array of spectrum object to velocity space
    spectrum.calculate_velocity()
    # find RA & DEC of target from target list
    spectrum.get_coords(TARGETS)
    # l & b are properties of the spectrum object dependent on RA & DEC and
    # are therefore automatically available
    # calculate velocity correction & add velocity correction to spectrum
    # object wavelength array (in velocity space)
    spectrum.lsr_correct_velocity()

    return spectrum


# --------------------------------------------------------------------------- #


def continuum_fit(spectrum, left=DEFAULTS["continuum_left"],
                  right=DEFAULTS["continuum_right"]):
    """This function performs all the calculations relating to the
    continuum, as well as calculating the normalized flux
    and error arrays"""

    # find indices in spectrum.velocity corresponding to left and right
    # continuum window boundaries
    left_indices = spectrum.find_indices(left)
    right_indices = spectrum.find_indices(right)
    indices = [left + right for left, right in
               zip(left_indices, right_indices)]
    # calculating means, signal to noise, pixel size, and S/N per resel
    continuum, signalnoise, pixels = spectrum.calculate_continuum(indices)

    # calculating fits to continuum windows and updating normalized flux and
    # error attributes
    spectrum, linear_fits = spectrum.calculate_fits(continuum)

    return spectrum


# --------------------------------------------------------------------------- #


def measure(args, spectrum):
    # find indices in velocity window (from -100 to 100, for example)
    window = [args.vel_min, args.vel_max]
    indices = spectrum.find_indices(window)
    # make truncated spectrum object to perform these measurements
    helper = Helper(spectrum, indices)
    # fix negative zeroes in error
    helper.error = helper.fix_negatives(helper.error)
    # set negative values in flux to percentage of the continuum array (will
    # mark these as saturation)
    helper.flux = helper.fix_negatives(helper.flux, fixes="negatives")
    # calculate apparent optical depth and error
    helper.calculate_aod()
    # calculate apparent column density and error
    helper.calculate_acd()
    # calculate equivalent width and error
    helper.calculate_ew()
    # evaluate significance of measurement, set detection flag
    helper.significance()

    # sets measurements done in helper object in spectrum object
    spectrum.set_measurements(helper)

    return spectrum


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    status = main()
    sys.exit(status)
