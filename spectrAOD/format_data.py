#!/usr/bin/env python

"""This module is for formatting input data into an object that stores primarily the wavelength and flux data, as well
as any relevant target information. This object is then used by measure_aod.py to perform the measurements."""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import os
import glob
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy import constants
from astropy.coordinates import SkyCoord
from operator import truediv

DATADIR = "/user/cmagness/fermi/data/15339/UVQSJ191928-295808/x1d/"

C = constants.c.to('km/s')  # km/s

N = 3.768e14  # proportionality constant -> (m_e * c)/(pi * e**2)

OUTDIR = "/user/cmagness/fermi/data/out/"  # this is unused at this time but will be eventually

# at some point need to check and see if outdir (and datadir, really too) exist. if not, create outdir
# add logger instead of print statements

# notes on class functionality:
# should have attributes for wavelength, flux at a minimum
# also would be great to have target name, l, b coordinates. or conversion from ra & dec if that information is
# available
# add functionality to build fits files as a method of class objects

# ----------------------------------------------------------------------------------------------------------------------


class Helper:
    """This class handles measuring AOD/ACD/EW operations"""

    def __init__(self, spectrum, indices):
        self.ions = spectrum.ions
        truncated_dict = {"velocity": [], "flux": [], "error": [], "continuum": [], "continuum_error": [], "delta": []}
        for idx, velocity in spectrum.velocity:
            truncated_dict["velocity"].append(velocity[indices[0]:indices[1]])
            truncated_dict["flux"].append(spectrum.flux[indices[0]:indices[1]])
            truncated_dict["error"].append(spectrum.norm_error[indices[0]:indices[1]])
            truncated_continuum = spectrum.continuum[idx][indices[0]:indices[1]]
            truncated_dict["continuum"].append(truncated_continuum)
            continuum_error = 0.05 * truncated_continuum  # 5% continuum fitting error
            truncated_dict["continuum_error"].append(continuum_error)
            truncated_dict["delta"].append(spectrum.delta[idx][indices[0]:indices[1]])
        self.velocity = truncated_dict["velocity"]
        self.flux = truncated_dict["flux"]
        self.error = truncated_dict["error"]
        self.continuum = truncated_dict["continuum"]
        self.cont_error = truncated_dict["continuum_error"]
        self.delta = truncated_dict["delta"]
        self.aod = None  # this will become a dictionary of values and errors
        self.acd = None
        self.ew = None
        self.sig = None
        self.detection = None

    def fix_negatives(self):
        # this method fixes any negative (and therefore unphysical) flux values
        for idx, subflux in enumerate(self.flux):
            if any(val < 0 for val in subflux):
                self.flux[idx] = [self.continuum[idx][subflux.index(val)] * 0.01 if val < 0 else val for val in subflux]

    def calculate_aod(self):
        # this method calculates the apparent optical depth and error
        aod = []
        for idx, subflux in enumerate(self.flux):
            aod_row = np.log(self.continuum[idx]/subflux)
            cont_err = map(truediv, self.cont_error[idx], self.continuum[idx])
            flux_err = map(truediv, self.error, subflux)  # FIX THIS
            total_err = np.sqrt([cont**2 for cont in cont_err] + [flux**2 for flux in flux_err])

            aod.append({"aod": [aod_row], "error": [total_err]})
        self.aod = aod

    @property
    def aod_val(self):
        # property to calculate single AOD sum for each ion measurement
        if self.aod:
            aod_val = []
            for idx, subdict in enumerate(self.aod):
                single_aod = np.sum([aod * delta for aod, delta in zip(subdict["aod"], self.delta[idx])])
                # each val in aod array is being multiplied by bin width
                single_err = np.sqrt(np.sum([(err * delta)**2 for err, delta in zip(subdict["aod"], self.delta[idx])]))

                aod_val.append({"aod": single_aod, "error": single_err})
        else:
            aod_val = None
        return aod_val

    def calculate_acd(self):
        # this method calculates the apparent column density and error
        acd = []
        for idx, subaod in enumerate(self.aod):
            acd_row = N/(self.ions["wv"][idx] * self.ions["f"][idx]) * subaod["aod"]
            acd_err = N/(self.ions["wv"][idx] * self.ions["f"][idx]) * subaod["error"]

            acd.append({"acd": [acd_row], "error": [acd_err]})
        self.acd = acd

    @property
    def acd_val(self):
        # property to calculate single ACD sum for each ion measurement
        if self.acd:
            acd_val = []
            for idx, subdict in enumerate(self.acd):
                linear_acd = np.sum([acd * delta for acd, delta in zip(subdict["acd"], self.delta[idx])])
                linear_err = np.sqrt(np.sum([(err * delta)**2 for err, delta in zip(subdict["error"], self.delta[idx])]))
                if linear_acd > 0:
                    log_acd = np.log10(linear_acd)
                else:
                    log_acd = 0.00
                log_err = np.log10(linear_acd + linear_err) - np.log10(linear_acd)

                acd_val.append(dict(linear_acd=linear_acd, linear_error=linear_err, log_acd=log_acd, log_error=log_err))
        else:
            acd_val = None
        return acd_val

    def calculate_ew(self):
        # this method calculates the equivalent width and error
        ew = []
        for idx, subflux in enumerate(self.flux):
            ew_row = (self.ions["wv"][idx]/C) * self.delta * (1.0 - map(truediv, subflux, self.continuum[idx]))
            cont_err = self.delta * subflux * (self.cont_error[idx] / self.continuum[idx] ** 2)  # FIX THIS AND CHECK
            flux_err = self.delta * map(truediv, self.error, subflux)  # FIX THIS
            total_err = (self.ions["wv"][idx]/C) * np.sqrt([cont ** 2 for cont in cont_err] + [flux ** 2 for flux in
                                                                                               flux_err])

            ew.append({"ew": [ew_row], "error": [total_err]})
        self.ew = ew

    @property
    def ew_val(self):
        # property to calculate single EW sum for each ion measurement
        if self.ew:
            ew_val = []
            for idx, subdict in enumerate(self.ew):
                single_ew = np.sum(subdict["ew"])
                single_err = np.sqrt(np.sum([err**2 for err in subdict["error"]]))

                ew_val.append({"ew": single_ew, "error": single_err})
        else:
            ew_val = None
        return ew_val

    def significance(self):
        # this method determines the significance of the measurement, checks for saturation & non detections
        sig = []
        for subdict in self.ew_val:
            sig.extend((subdict["ew"] >= 3.0 * subdict["error"]))
        self.sig = sig

        # check that this is what we want the logic to be for these detection moments
        detection = []
        for idx, subflux in enumerate(self.flux):
            minimum = min(subflux/self.continuum[idx])
            cutoff = 0.1
            if minimum > cutoff:
                if self.sig[idx]:
                    detection.extend("D")
                else:
                    detection.extend("ND")
            else:  # minimum <= cutoff
                detection.extend("S")
        self.detection = detection

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
        self.sn_avg = None
        self.sn_res = None
        self.delta = None
        self.ions = None
        self.doublet = False
        self.velocity = None
        self.ra = None
        self.dec = None
        self._skycoords = None

        # things that will be retrieved from helper class
        self.aod = None
        self.acd = None
        self.ew = None
        self.aod_val = None
        self.acd_val = None
        self.ew_val = None
        self.sig = None
        self.detection = None
        # ADD A FLAG FOR INSIDE VS OUTSIDE FERMI BUBBLES

    @property
    def l(self):
        return self._skycoords.galactic.l.value

    @property
    def b(self):
        return self._skycoords.galactic.b.value

    def get_ions(self, ion, file="mini_ions.csv"):
        # this method retrieves the ion of interest's wavelength from the ions file and creates a dictionary. Accounts
        # for doublets as well.
        df_ions = pd.read_csv(file, delimiter=" ", header=None)
        df_masked = df_ions[df_ions[0] == ion]
        # DO I EVEN NEED THE DOUBLET FLAG
        doublet = False
        if len(df_masked) > 1:
            doublet = True
        ions = {"ion": [], "wv": [], "f": []}
        for index, row in df_masked.iterrows():
            ion = row[0]
            wv = float(row[1])  # truncate this to 4 digits
            f_value = float(row[2])
            ions["ion"] += [ion]
            ions["wv"] += [wv]
            ions["f"] += [f_value]

        self.doublet = doublet
        self.ions = ions

    def calculate_velocity(self):
        # this method calculates the wavelength in velocity space wrt to the ion of interest and stores it
        for wv in self.ions["wv"]:
            z_array = (self.wave - wv) / wv
            vel_array = C * z_array
            self.velocity.append(vel_array)
        # eventually will need to keep track of ion associated with the velocity array in a better way maybe?
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

    def find_indices(self, window):
        # this method finds the indices for a velocity window for each velocity array
        indices = []
        for idx in np.arange(len(self.velocity) + 1):
            index_row = []
            for val in window:
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

            upper = self.velocity[idx][2:]
            lower = self.velocity[idx][:-2]
            delta = (upper - lower)/2.0
            dv_row = [delta[0], delta, delta[-1]]
            dv.append([dv_row])

        self.sn_avg = [sn_row[-1] for sn_row in signalnoise]
        self.delta = dv

        pixels = []
        for idx, dv_row in enumerate(dv):
            index_row = indices[idx]
            pixsize = np.mean(dv_row[index_row[0]:index_row[1]])
            sn_res = signalnoise[idx][2] * np.sqrt(
                C / (16000.0 * pixsize))  # signal to noise per resolution element
            pixels.append([pixsize, sn_res])

        self.sn_res = [sn_row[-1] for sn_row in pixels]

        # SHOULD PROBABLY CLEAN UP ALL THE UNUSED STUFF IN THIS METHOD

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

    def set_measurements(self, helper):
        # this method will set the measurements as calculated by the helper object as attributes in this object
        self.aod = helper.aod
        self.aod_val = helper.aod_val
        self.acd = helper.acd
        self.acd_val = helper.acd_val
        self.ew = helper.ew
        self.ew_val = helper.ew_val
        self.sig = helper.sig
        self.detection = helper.detection

    def generate_table(self, vel_min, vel_max):
        # this method should generate the table with all the measurements
        cols = ["TARGET", "RA", "DEC", "ION", "VEL MIN", "VEL MAX", "AOD", "AOD ERR", "ACD (LIN)", "ACD ERR (LIN)",
                "ACD (LOG)", "ACD ERR (LOG)", "EW", "EW ERR", "SN AVG", "SN/RESEL", "SIG", "DETECTION"]
        df = pd.DataFrame(columns=cols)
        for idx, row in enumerate(self.velocity):
            df.append({
                "TARGET": self.target,
                "RA": self.ra,
                "DEC": self.dec,
                "ION": self.ions[idx],
                "VEL MIN": vel_min,
                "VEL MAX": vel_max,
                "AOD": self.aod_val[idx]["aod"],  # ok really the _val properties COULD be properties of this class
                "AOD ERR": self.aod_val[idx]["error"],  # instead of properties of the helper class
                "ACD (LIN)": self.acd_val[idx]["linear_acd"],  # which would eliminate the need to pass them between
                "ACD ERR (LIN)": self.acd_val[idx]["linear_error"],  # consider this in refactoring
                "ACD (LOG)": self.acd_val[idx]["log_acd"],
                "ACD ERR (LOG)": self.acd_val[idx]["log_error"]
                "EW": self.ew_val[idx]["ew"],
                "EW ERR": self.ew_val[idx]["error"],
                "SN AVG": self.sn_avg[idx],
                "SN/RESEL": self.sn_res[idx],
                "SIG": self.sig[idx],
                "DETECTION": self.detection[idx]
            })
        df.to_csv(os.path.join(OUTDIR + "measurements.csv"))

    def to_fits(self):
        # this method should generate a fits file
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
