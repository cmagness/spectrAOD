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


def continuum_fit(args, spectrum, left=[-450, -300], right=[300, 450]):

    # find indices in spectrum.velocity corresponding to left and right continuum boundaries
    indices = []
    for idx in np.arange(len(spectrum.velocity) + 1):
        index_row = []
        for val in left[0], left[1], right[0], right[1]:
            index_row.append(list(spectrum.velocity[idx]).index(val))
        indices.append(index_row)

    # continuum will have flux, velocity, and noise
    continuum = []  # this is a list holding the continuum measurements for each spectrum.velocity
    for idx, row in enumerate(indices):
        left_min, left_max, right_min, right_max = row
        # substitute below
        flux_l = np.mean(spectrum.flux[row[0]:row[1]])  # mean of flux from in range of indices determined above
        flux_r = np.mean(spectrum.flux[row[2]:row[3]])
        vel_l = np.mean(spectrum.velocity[idx][row[0]:row[1]])  # mean of velocity array corresponding to row in range
        vel_r = np.mean(spectrum.velocity[idx][row[2]:row[3]])  # of indices determined above
        noise_l = np.std(spectrum.flux[row[0]:row[1]])   # std dev of flux in range as determined above
        noise_r = np.std(spectrum.flux[row[2]:row[3]])
        continuum_row = {"flux": [flux_l, flux_r], "velocity": [vel_l, vel_r], "noise": [noise_l, noise_r]}
        # ^ this can be stored as left and right instead of by flux, velocity, and noise
        continuum.append(continuum_row)

    # do we want any of this information ^^ in the final table? if so we just need to add the continuum list to the
    # Spectrum object
    # i don't think so based on discussions

    signalnoise = []
    for idx, row in enumerate(indices):
        continuum_row = continuum[idx]
        sn_l = continuum_row["flux"][0]/continuum_row["noise"][0]
        sn_r = continuum_row["flux"][1]/continuum_row["noise"][1]
        sn_avg = (sn_l + sn_r)/2.0
        sn_row = [sn_l, sn_r, sn_avg]
        signalnoise.append(sn_row)

    # what are we doing with this signal to noise

    fit = []
    for idx, cdict in enumerate(continuum):
        slope = (cdict["flux"][1] - cdict["flux"][0])/(cdict["velocity"][1] - cdict["velocity"][0])  # right - left f/v
        yint = cdict["flux"][0] - (slope * cdict["velocity"][0])
        continuum_array = (slope * spectrum.velocity[idx]) + yint
        normalized_flux = spectrum.flux/continuum_array
        normalized_error = spectrum.error/continuum_array  # DO WE NEED THIS BC IF SO I NEED TO GET ERROR FROM THE DATA
        fit_row = [slope, yint, continuum_array, normalized_flux, normalized_error]
        fit.append(fit_row)

    # what are we doing with all this fit stuff ^^

    dv = []
    for velocity in spectrum.velocity:
        upper = velocity[1:]
        lower = velocity[:-1]
        deriv = upper - lower
        dv_row = [deriv[0], deriv]
        dv.append(dv_row)

    # what are we doing with this ?? ^^

    pixels = []
    for idx, dv_row in enumerate(dv):
        index_row = indices[idx]
        pixsize = np.mean(dv_row[index_row[0]:index_row[1]])
        sn_res = signalnoise[idx][2] * np.sqrt(2.998e5/(16000.0 * pixsize))  # signal to noise per resolution element
        pixels.append([pixsize, sn_res])

    # also what are we doing with THIS ?? ^^

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
    # does this need a separate function?
    # generate table
    spectrum.generate_table()
    return 0


# ----------------------------------------------------------------------------------------------------------------------


if __name__ == "__main__":
    status = main()
    sys.exit(status)
