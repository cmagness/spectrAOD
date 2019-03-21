#!/usr/bin/env python

"""This module is for measuring the apparent optical depth of absorption lines in COS Fermi Bubble spectra. The method
of measurement is outlined in Savage & Sembach 1991. Eventually, this will be extended to measure apparent optical depth
of other spectra as well."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import argparse
import sys

from .format_data import *

DATADIR = "/user/cmagness/fermi/data/15339/UVQSJ191928-295808/x1d/"
OUTDIR = "/user/cmagness/fermi/data/out/"

# ----------------------------------------------------------------------------------------------------------------------


def parse():
    parser = argparse.ArgumentParser(description="spectrAOD")

    parser.add_argument("instrument", type=str, help="observational instrument data is taken on")
    parser.add_argument("filetype", type=str, help="file type of data")
    parser.add_argument("ion", type=str, help="absorption line feature of interest")
    parser.add_argument("velocity", type=int, help="velocity window of interest around ion, in km")
    parser.add_argument("--grating", type=str, help="optional grating to choose from")

    # grating evaluation needs some work because it will have to be checked that the ion of interest falls in wavelength
    # range of the grating

    args = parser.parse_args()

    return args

# ----------------------------------------------------------------------------------------------------------------------


def lsr_correct(args):
    # helio to LSR correction in km / s:
    # longr = long[i] *!pi / 180. & latr = lat[i] *!pi / 180.
    # vcorrd = 9. * cos(longr) * cos(latr) + 12. * sin(longr) * cos(latr) + 7. * sin(latr)

    # file with lat & long of targets is at /user/afox/research/gc/sl.cos.dat which i do not have access to but
    # will need. will need to cross check that with target list?
    # how does this correction get applied to the Spectrum object?

    spectrum = format_data(DATADIR, args.instrument, args.filetype, args.grating)

    # this is just some placeholder code
    lsr = 1
    corrected_spectrum = spectrum + lsr

    return corrected_spectrum


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
    args = parse()
    # LSR correction
    spectrum = lsr_correct(args)
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
