#!/usr/bin/env python

"""This module is for measuring the apparent optical depth of absorption lines in COS Fermi Bubble spectra. The method
of measurement is outlined in Savage & Sembach 1991. Eventually, this will be extended to measure apparent optical depth
of other spectra as well."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import sys
import argparse

from .format_data import *

C = 2.99792458e5  # km/s

DATADIR = "/user/cmagness/fermi/data/15339/UVQSJ191928-295808/x1d/"
OUTDIR = "/user/cmagness/fermi/data/out/"

TARGETS = "/user/cmagness/fermi/code/agn_target_list.csv"

# ----------------------------------------------------------------------------------------------------------------------


def parse():
    """This function is for parsing the inputs of interest."""

    parser = argparse.ArgumentParser(description="spectrAOD")

    parser.add_argument("instrument", type=str, help="observational instrument data is taken on")
    parser.add_argument("filetype", type=str, help="file type of data")
    parser.add_argument("ion", type=str, help="absorption line feature of interest")
    parser.add_argument("velocity", type=int, help="velocity window of interest around ion, in km")
    parser.add_argument("--grating", type=str, help="optional grating to choose from")
    # need to set default grating to something for creating Spectrum object

    # grating evaluation needs some work because it will have to be checked that the ion of interest falls in wavelength
    # range of the grating

    args = parser.parse_args()

    spectrum = format_data(DATADIR, args.instrument, args.filetype, args.grating)

    return args, spectrum

# ----------------------------------------------------------------------------------------------------------------------


def get_ion(ion, file="ions.csv"):
    return ion_wv
    pass

# ----------------------------------------------------------------------------------------------------------------------


def lsr_correct(args, spectrum):
    """This function performs the lsr correction"""

    # find ion wavelength from ions.csv
    ion_wv = get_ion(args.ion)
    # transform wavelength array of Spectrum object to velocity space
    spectrum.calculate_velocity(ion_wv)
    # find RA & DEC of target from target list
    spectrum.get_coords(TARGETS)
    # transform RA & DEC to L & B
    spectrum.to_galactic()
    # calculate velocity correction & add velocity correction to Spectrum object wavelength array (in velocity space)
    spectrum.lsr_correct_velocity()

    # to check: do i need to reassign the spectrum to itself since i've updated attributes?

    return spectrum


# ----------------------------------------------------------------------------------------------------------------------


def continuum_fit(args, spectrum):
    pass


# ----------------------------------------------------------------------------------------------------------------------


def measure_aod():

    pass


# ----------------------------------------------------------------------------------------------------------------------


def main():
    # steps
    # argparse
    args, spectrum = parse()
    # LSR correction
    spectrum = lsr_correct(args, spectrum)
    # continuum fit
    continuum_fit(spectrum)
    # measure aod
    measure_aod()
    # store measurements
    return 0


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    status = main()
    sys.exit(status)
