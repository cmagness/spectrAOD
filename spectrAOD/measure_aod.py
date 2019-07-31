#!/usr/bin/env python

"""This module is for measuring the apparent optical depth of absorption lines
in COS Fermi Bubble spectra. The method of measurement is outlined in Savage &
Sembach 1991. Eventually, this will be extended to measure apparent optical
depth of other spectra as well."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import argparse
import sys

from format_data import *

# from astropy import constants

# C = constants.c.to('km/s')  # km/s

DATADIR = "/user/cmagness/fermi/data/15339/UVQSJ191928-295808/x1d/"
OUTDIR = "/user/cmagness/fermi/data/out/"

TARGETS = "/user/cmagness/fermi/code/agn_target_list.csv"


# -----------------------------------------------------------------------------


def parse():
    """This function is for parsing the inputs of interest."""

    parser = argparse.ArgumentParser(description="spectrAOD")

    parser.add_argument("--instrument", default="COS", type=str,
                        help="observational instrument data is taken on")
    parser.add_argument("--filetype", default="X1DSUM", type=str,
                        help="file type of data")
    parser.add_argument("--ion", default="SiIV", type=str,
                        help="absorption line feature of interest")
    parser.add_argument("--vel_min", default="-100", type=int,
                        help="velocity minimum in window of interest around "
                             "ion, in km/s")
    parser.add_argument("--vel_max", default="100", type=int,
                        help="velocity maximum in window of interest around "
                             "ion, in km/s")
    parser.add_argument("--grating", default="G130M", type=str,
                        help="optional grating to choose from")
    # need to set up better debugging system than defaults. debugging flag?
    # need to set default grating to something for creating Spectrum object
    # NEED TO DO SOMETHING ABOUT THE FACT THAT COULD HAVE ION LIST IN THE
    # FUTURE

    # grating evaluation needs some work because it will have to be checked
    # that the ion of interest falls in wavelength
    # range of the grating

    args = parser.parse_args()
    # a = vars(args)
    if args.grating:
        spectrum = collect(DATADIR, args.instrument, args.filetype,
                           args.grating)
    else:
        spectrum = collect(DATADIR, args.instrument, args.filetype)

    # switch the arguments to unpacking with vars() and switch collect to
    # accepting **kw args

    return args, spectrum


# -----------------------------------------------------------------------------


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


# -----------------------------------------------------------------------------


def continuum_fit(spectrum, left=[-450, -300], right=[300, 450]):
    """This function performs all the calculations relating to the
    continuum, as well as calculating the normalized flux
    and error arrays"""

    # find indices in spectrum.velocity corresponding to left and right
    # continuum window boundaries
    left_indices = spectrum.find_indices(left)
    right_indices = spectrum.find_indices(right)
    indices = [left + right for left, right in zip(left_indices,
                                                   right_indices)]
    # calculating means, signal to noise, pixel size, and S/N per resel
    continuum, signalnoise, pixels = spectrum.calculate_continuum(indices)
    # WE WANT TO ADD THE SN_AVG AS AN ATTRIBUTE FOR THE TABLE AT SOME POINT
    # PROBABLY
    # calculating fits to continuum windows and updating normalized flux and
    # error attributes

    # note: pixels is empty at this point after the test run
    # this needs investigating
    spectrum, linear_fits = spectrum.calculate_fits(continuum)

    return spectrum


# -----------------------------------------------------------------------------


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


def main():
    # argparse
    args, spectrum = parse()
    # LSR correction
    spectrum = lsr_correct(args, spectrum)
    # continuum fit
    spectrum = continuum_fit(spectrum)
    # measure aod/acd/ew
    # set measurements back in spectrum object from helper object
    spectrum = measure(args, spectrum)
    # generate table
    spectrum.generate_table(args.vel_min, args.vel_max)
    return 0


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    status = main()
    sys.exit(status)
