#!/usr/bin/env python

"""This module is for running the measurement process multiple times.
Requires the use of batch_table.csv"""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import sys
import glob
import logging

import pandas as pd
from tqdm import tqdm
from astropy import constants

from . import SETTINGS
from .measure_aod import *
from .format_data import build_spectrum
from .visualization import Visualizer

INPUTS = SETTINGS["inputs"]
DATADIR = INPUTS["datadir"]
PARAMETERS = SETTINGS["parameters"]
DEFAULTS = SETTINGS["defaults"]
LOGGER = logging.getLogger(__name__)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
LOGGER.addHandler(console)

C = constants.c.to('km/s').value  # km/s

N = 3.768e14  # proportionality constant -> (m_e * c)/(pi * e**2)


# at some point need to check and see if outdir (and datadir, really too)
# exist. if not, create outdir


# --------------------------------------------------------------------------- #


def main():
    # look at the table for info-turn into df, return df
    LOGGER.info("Entering batch mode...")
    if DEFAULTS["batch_table"]:
        batch_dataframe = read_table(DEFAULTS["batch_table"])
    else:
        LOGGER.error("Please provide a path to your batch table in your "
                     "settings file if you would like to perform multiple "
                     "measurements. Exiting...")
        sys.exit()

    # use rootnames to collect the files that need to be run on
    batch_dataframe = collect_files(batch_dataframe)

    # for each file, essentially run a slightly modified version of
    # measure_aod.main
    batch_run(batch_dataframe)
    return 0


def read_table(filename):
    # why did i make this its own function lol
    return pd.read_csv(filename)


def collect_files(dataframe):
    # get the column of the basenames as a list or w/e
    rootnames = dataframe["ROOTNAME"]

    # get the files in the data directory
    all_files_in_dir = glob.glob(DATADIR + "*")

    # for each file in the column, check to see if a file matches
    batch_files = []
    for rootname in rootnames:
        for filename in all_files_in_dir:
            if rootname in filename:
                batch_files += [filename]
                break
        else:  # only gets here if it doesn't break the filename statement,
            # i.e. it didn't find a match
            LOGGER.warning("No file was found matching the rootname: {}. "
                           "Continuing...".format(rootname))

    if batch_files:
        LOGGER.info("Found {} files to measure. This might take a while."
                    .format(len(batch_files)))
        dataframe["FILENAME"] = batch_files
    else:
        LOGGER.warning("Found no files to measure. Exiting...")
        sys.exit()

    return dataframe


def batch_run(dataframe):
    # for each file in the list, do:
    for index, file_row in tqdm(dataframe.iterrows()):
        # collect the arguments
        args = {"datadir": file_row["FILENAME"],
                "ins": file_row["INSTRUMENT"].upper(),
                "file": file_row["FILETYPE"].upper(),
                "grating": file_row["GRATING"].upper(),
                "redshift": file_row["REDSHIFT"]}

        # build a spectrum object
        spectrum = build_spectrum(**args)
        spectrum.target = file_row["TARGET"]  # this is a hack to get around
        # stuff in build_spectrum
        LOGGER.info("Spectrum object successfully built.")

        # pass that on and do everything exactly the same
        # set up visualizer
        visualizer = Visualizer()
        visualizer.set_target(spectrum.target)
        visualizer.set_raw_flux(spectrum.flux)

        # LSR correction
        # need to add ion in to the args here
        args["ion"] = file_row["ION"]
        spectrum = lsr_correct(args, spectrum)
        LOGGER.info("Spectrum LSR corrected.")
        visualizer.set_raw_velocity(spectrum.raw_velocity[0])
        visualizer.set_lsr_velocity(spectrum.velocity[0])

        # continuum fit
        spectrum, left_indices, right_indices = continuum_fit(spectrum)
        LOGGER.info("Continuum fit calculated.")
        visualizer.set_contadjspec(spectrum)
        visualizer.set_indices(left_indices, right_indices)

        # measure aod/acd/ew
        # set measurements back in spectrum object from helper object
        # need to add the vel_min & max here
        args["vel_min"] = file_row["VEL_MIN"]
        args["vel_max"] = file_row["VEL_MAX"]
        spectrum, helper = measure(args, spectrum)
        visualizer.set_helper(helper)

        # generate table
        spectrum.generate_table(args["vel_min"], args["vel_max"])
        visualizer.plot()
        LOGGER.info("Finished measurements for {}"
                    .format(file_row["FILENAME"]))

    # finish entire list
    LOGGER.info("spectrAOD complete.")


# --------------------------------------------------------------------------- #


if __name__ == "__main__":
    status = main()
    sys.exit(status)
