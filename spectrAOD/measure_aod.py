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


def continuum_fit(args, spectrum):
    # linear fit between continuum velocity ranges
    # how should we choose these ranges?
    # find index of velocity ranges
    # continuum will have flux, velocity, and noise
    # flux = mean of flux between index window
    # velocity = mean of velocity between index window
    # noise = standard deviation of the flux between index window
    # S/N = continuum/noise for both continuum ranges
    # avg S/N

    # slope=(cont2-cont1)/(vc2-vc1) & yint=cont1-slope*vc1
    # continuum array = slope*v + yint & fnorm=flux/carr & enorm=error/carr
    # dv=v & for k=1, n_elements(v)-1 do dv[k]=v[k]-v[k-1] & dv[0]=dv[1]
    # pixsize=mean(dv[i1:i2]) & stonres=stona*sqrt(2.998d5/(16000.*pixsize)) ; S/N per resolution element
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
