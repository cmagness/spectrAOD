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
    # NEED TO DO SOMETHING ABOUT THE FACT THAT COULD HAVE ION LIST IN THE FUTURE

    # grating evaluation needs some work because it will have to be checked that the ion of interest falls in wavelength
    # range of the grating

    args = parser.parse_args()

    spectrum = format_data(DATADIR, args.instrument, args.filetype, args.grating)

    return args, spectrum

# ----------------------------------------------------------------------------------------------------------------------


def get_ion(ion, spectrum, file="mini_ions.csv"):
    """This function retrieves the ion of interest's wavelength from the ions file and creates a dictionary. Accounts for
    doublets as well."""

    df_ions = pd.read_csv(file, delimiter=" ", header=None)
    df_masked = df_ions[df_ions[0] == ion]
    # DO I EVEN NEED THE DOUBLET FLAG
    doublet = False
    if len(df_masked) > 1:
        doublet = True
    ion_wv = {}
    for index, row in df_masked.iterrows():
        row_ion = row[0]
        wavelength = float(row[1])  # truncate this to 4 digits
        ion_wv[wavelength] = row_ion

    spectrum.set_doublet(doublet)

    return ion_wv

# ----------------------------------------------------------------------------------------------------------------------


def lsr_correct(args, spectrum):
    """This function performs the lsr correction"""

    # find ion wavelength from ions.csv
    ion_wv = get_ion(args.ion, spectrum)
    # transform wavelength array of Spectrum object to velocity space
    spectrum.calculate_velocity(ion_wv)
    # find RA & DEC of target from target list
    spectrum.get_coords(TARGETS)
    # l & b are properties of the Spectrum object dependent on RA & DEC and are therefore automatically available
    # calculate velocity correction & add velocity correction to Spectrum object wavelength array (in velocity space)
    spectrum.lsr_correct_velocity()

    return spectrum


# ----------------------------------------------------------------------------------------------------------------------


def continuum_fit(spectrum, left=[-450, -300], right=[300, 450]):
    """This function performs all the calculations relating to the continuum, as well as calculating the normalized flux
    and error arrays"""

    # find indices in spectrum.velocity corresponding to left and right continuum window boundaries
    indices = find_indices(spectrum, left, right)
    # calculating means, signal to noise, pixel size, and S/N per resel
    continuum, signalnoise, pixels = calculate_continuum(spectrum, indices)
    # calculating fits to continuum windows and updating normalized flux and error attributes
    spectrum, fits = calculate_fits(spectrum, continuum)

    return spectrum

# ----------------------------------------------------------------------------------------------------------------------


def measure_aod():

    # find indices in velocity window (from -100 to 100)
    # truncate the full arrays into sub arrays around each window
    # don't calculate pixel size again
    # set negative values in flux to percentage of the continuum array (will mark these as saturation)
    # aod is the natural log of the continuum / flux which is 1 / normalized flux
    # calculate errors on aod  (velocity error & flux error), then add in quadrature
    # calculate total aod by summing every element in array
    # same set of calculations for acd
    # f val and lambda come from the ions.csv, wavelength and f value

    # EW
    pass


# ----------------------------------------------------------------------------------------------------------------------


def main():
    # steps
    # argparse
    args, spectrum = parse()
    # LSR correction
    spectrum = lsr_correct(args, spectrum)
    # continuum fit
    spectrum = continuum_fit(spectrum)
    # measure aod
    measure_aod()
    # store measurements
    # does this need a separate function?
    # also perhaps a separate function for measuring the equivalent width
    # generate table
    spectrum.generate_table()
    return 0


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    status = main()
    sys.exit(status)
