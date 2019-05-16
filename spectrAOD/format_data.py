#!/usr/bin/env python

"""This module is for formatting input data into an object that stores primarily the wavelength and flux data, as well
as any relevant target information. This object is then used by measure_aod.py to perform the measurements."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import glob
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

DATADIR = "/user/cmagness/fermi/data/15339/UVQSJ191928-295808/x1d/"

C = 2.99792458e5  # km/s


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
    """This is the base class for the Spectrum objects that get used to store data for this package."""

    # NEED TO ADD SOME ERROR HANDLING IF METHODS ARE CALLED BEFORE ATTRIBUTES USED ARE DEFINED

    def __init__(self, target, wave, flux, error):
        self.target = target
        self.wave = wave
        self.flux = flux
        self.error = error
        self.norm_flux = None
        self.norm_error = None
        self.continuum = None
        self.doublet = False
        self.velocity = None
        self.ra = None
        self.dec = None
        self._skycoords = None
        # ADD A FLAG FOR INSIDE VS OUTSIDE FERMI BUBBLES

    @property
    def l(self):
        return self._skycoords.galactic.l.value

    @property
    def b(self):
        return self._skycoords.galactic.b.value

    def set_doublet(self, doublet):
        self.doublet = doublet

    def calculate_velocity(self, ion_wv: dict):
        # this method calculates the wavelength in velocity space wrt to the ion of interest and stores it
        for wv in ion_wv:
            z_array = (self.wave - wv) / wv
            vel_array = C * z_array
            self.velocity.append(vel_array)
        # eventually will need to keep track of ion associated with the velocity array
        # for now it is fine because they are associated with the same ion in the doublet case

    def get_coords(self, target_list):
        # this method should get the RA & DEC from the coordinates list
        df_targets = pd.read_csv(target_list, index_col=0)
        # need to add some error handling for if target name is not found in list
        mask = df_targets["Target"] == self.target
        self.ra = df_targets.loc[mask]["RA"]
        self.dec = df_targets.loc[mask]["DEC"]
        self._skycoords = SkyCoord(ra=self.ra * u.degree, dec=self.dec * u.degree)

    def lsr_correct_velocity(self):
        # this method should apply the lsr correction to the heliocentric coordinates
        l_radians = self.l * np.pi / 180.0
        b_radians = self.b * np.pi / 180.0
        vel_corr = (9.0 * np.cos(l_radians) * np.cos(b_radians)) + (12.0 * np.sin(l_radians) * np.cos(b_radians)) + \
                   (7.0 * np.sin(b_radians))
        for idx in np.arange(len(self.velocity)):
            self.velocity[idx] = self.velocity[idx] + vel_corr

    def find_indices(self, left, right):
        # this method finds the indices for the continuum windows
        indices = []
        for idx in np.arange(len(self.velocity) + 1):
            index_row = []
            for val in left[0], left[1], right[0], right[1]:
                index_row.append(list(self.velocity[idx]).index(val))
            indices.append(index_row)
        return indices

    def calculate_continuum(self, indices):
        # this method calculates the continuum, signal to noise, and pixel measurements
        continuum = []  # this is a list holding the continuum measurements for each spectrum.velocity
        signalnoise = []
        dv = []
        for idx, row in enumerate(indices):
            left_min, left_max, right_min, right_max = row
            flux_l = np.mean(self.flux[left_min:left_max])  # mean of flux from in range of indices determined above
            flux_r = np.mean(self.flux[right_min:right_max])
            vel_l = np.mean(
                self.velocity[idx][left_min:left_max])  # mean of velocity array corresp. to row in range
            vel_r = np.mean(self.velocity[idx][right_min:right_max])  # of indices determined above
            noise_l = np.std(self.flux[left_min:left_max])  # std dev of flux in range as determined above
            noise_r = np.std(self.flux[right_min:right_max])
            continuum_row = {"flux": [flux_l, flux_r], "velocity": [vel_l, vel_r], "noise": [noise_l, noise_r]}
            # ^ this can be stored as left and right instead of by flux, velocity, and noise
            continuum.append(continuum_row)

            sn_l = continuum_row["flux"][0] / continuum_row["noise"][0]
            sn_r = continuum_row["flux"][1] / continuum_row["noise"][1]
            sn_avg = (sn_l + sn_r) / 2.0
            sn_row = [sn_l, sn_r, sn_avg]
            signalnoise.append(sn_row)

            upper = self.velocity[idx][1:]
            lower = self.velocity[idx][:-1]
            delta = upper - lower
            dv_row = [delta[0], delta]
            dv.append(dv_row)

        pixels = []
        for idx, dv_row in enumerate(dv):
            index_row = indices[idx]
            pixsize = np.mean(dv_row[index_row[0]:index_row[1]])
            sn_res = signalnoise[idx][2] * np.sqrt(
                2.998e5 / (16000.0 * pixsize))  # signal to noise per resolution element
            pixels.append([pixsize, sn_res])

        return continuum, signalnoise, pixels  # figure out what to do with this later

    def calculate_fits(self, continuum):
        # this method calculates the linear fits to the continuum windows and calculates the normalized arrays
        linear_fits = []
        lists = []
        for idx, cdict in enumerate(continuum):
            slope = (cdict["flux"][1] - cdict["flux"][0]) / (
                    cdict["velocity"][1] - cdict["velocity"][0])  # right - left f/v
            yint = cdict["flux"][0] - (slope * cdict["velocity"][0])
            continuum_array = (slope * self.velocity[idx]) + yint
            normalized_flux = self.flux / continuum_array
            normalized_error = self.error / continuum_array
            fit_row = [slope, yint]
            list_row = [continuum_array, normalized_flux, normalized_error]
            linear_fits.append(fit_row)
            lists.append(list_row)
        self.continuum = lists[0]  # list of calculated continuums
        self.norm_flux = lists[1]  # list of normalized fluxes
        self.norm_error = lists[2]  # list of normalized errors

        return self, fits

    def to_fits(self):
        # this method should generate a fits file
        pass

    def generate_table(self):
        # this method should generate the table with all the measurements
        pass


# ----------------------------------------------------------------------------------------------------------------------


class X1DSpectrum(BaseSpectrum):
    """This class inherits the base Spectrum class for specifically x1dsum files and will have methods for holding
    other x1d specific information."""

    def __init__(self, *args):
        super().__init__(*args)

    def x1d_specs(self):
        # this method could hold other x1d specific stuff, like header information of interest
        pass

    # we're gonna wanna add other things we can capture with x1d spectra here later. which will also require updates to
    # that part of the if statement in format data


# ----------------------------------------------------------------------------------------------------------------------

# add classes for other data types

# ----------------------------------------------------------------------------------------------------------------------


def format_data(datadir=DATADIR, ins="COS", file="X1DSUM", grating="G130M"):
    # this function should build and return the appropriate Spectrum object for measure_aod.py
    # default inputs for instrument and file type are COS and X1D at the moment, for testing

    spectrum = []

    # this might (?) need handling for other inputs other than just COS as well. at some point.
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
                        wave = target_data["WAVELENGTH"]
                        flux = target_data["FLUX"]
                        error = target_data["ERROR"]
                        spectrum = X1DSpectrum(target, wave, flux, error)
        else:
            print("Other file types are not yet supported at this time.")
    else:
        print("This Instrument is not yet supported at this time.")

    # this needs to be addressed with proper error handling so that it is not referenced before assignment
    if not spectrum:
        spectrum = []
        print("Spectrum object not built.")

    return spectrum
