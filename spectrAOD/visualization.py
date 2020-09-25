#!/usr/bin/env python

"""Module Docstring"""

__author__ = "Camellia Magness"
__email__ = "cmagness@stsci.edu"

import os
import sys
import glob
import logging

import numpy as np
from matplotlib import rc
from astropy import constants
import matplotlib.pyplot as plt
from astropy.convolution import Box1DKernel, convolve

from . import SETTINGS

INPUTS = SETTINGS["inputs"]
DATADIR = INPUTS["datadir"]
OUTDIR = INPUTS["outdir"]
PARAMETERS = SETTINGS["parameters"]
DEFAULTS = SETTINGS["defaults"]
LOGGER = logging.getLogger(__name__)
console = logging.StreamHandler()
console.setLevel(logging.INFO)
LOGGER.addHandler(console)

C = constants.c.to('km/s').value  # km/s

N = 3.768e14  # proportionality constant -> (m_e * c)/(pi * e**2)


# --------------------------------------------------------------------------- #


class Visualizer:

    def __init__(self):
        self._target = None
        self._raw_flux = None
        self._raw_velocity = None
        self._lsr_velocity = None
        self._contadjspec = None  # this is the continuum adjusted spectrum
        self._left_indices = None
        self._right_indices = None
        self._helper = None

    def set_target(self, value):
        self._target = value

    def set_raw_flux(self, value):
        self._raw_flux = value

    def set_raw_velocity(self, value):
        self._raw_velocity = value

    def set_lsr_velocity(self, value):
        self._lsr_velocity = value

    def set_contadjspec(self, value):
        self._contadjspec = value

    def set_indices(self, left_value, right_value):
        self._left_indices = left_value
        self._right_indices = right_value

    def set_helper(self, value):
        self._helper = value

    def plot(self, box=5, show=DEFAULTS["show_plot"]):
        target = self._target
        raw_flux = self._raw_flux
        # at this point a list of arrays
        raw_velocities = self._raw_velocity
        # also a list of arrays
        lsr_velocities = self._lsr_velocity
        contadjspec = self._contadjspec
        helper = self._helper
        cont_left = DEFAULTS["continuum_left"]
        cont_right = DEFAULTS["continuum_right"]
        left = self._left_indices
        right = self._right_indices
        vel_min = PARAMETERS["vel_min"]
        vel_max = PARAMETERS["vel_max"]

        # you can't just do list membership, i.e. if None in <list>,
        # here because there are arrays in the list and that somehow messes
        # this up
        if any(elem is None for elem in [raw_flux, raw_velocities,
                                         lsr_velocities, contadjspec,
                                         helper]):
            LOGGER.error("NoneType in list: \n"
                         "Raw Flux: {} \n"
                         "Raw Velocity: {} \n"
                         "LSR Corrected Velocity: {} \n"
                         "Continuum Adjusted Spectrum: {} \n"
                         "Helper Object: {}".format(type(raw_flux),
                                                    type(raw_velocities),
                                                    type(lsr_velocities),
                                                    type(contadjspec),
                                                    type(helper)))
            LOGGER.error("Missing elements necessary to "
                         "plot. Exiting...")
            sys.exit()

        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')
        plt.rc('axes', axisbelow=True)

        for idx, ion in enumerate(helper.ions["ion"]):
            wv = helper.ions["wv"][idx]
            # need to extract the velocities that correspond to the correct ion
            raw_velocity = raw_velocities[idx]
            lsr_velocity = lsr_velocities[idx]
            fig, ax = plt.subplots(2, 1, figsize=(15, 10), gridspec_kw={
                'height_ratios': [2, 1]
                })
            fig.suptitle("Visualization for {}: {}-{}"
                         .format(target, ion, int(wv)),
                         fontsize=20)
            ax[0].set_xlim(-500, 500)
            ax[1].set_xlim(-500, 500)

            max_flux = max(raw_flux[left[idx][0]:right[idx][1]])
            ax[0].set_ylim(0, 2 * max_flux)
            # scaling the axis based on the maximum value from lower left
            # continuum bound to upper right continuum bound

            kernel = Box1DKernel(box)
            boxcar_flux = convolve(raw_flux, kernel)

            ax[0].plot(raw_velocity, boxcar_flux, color="cornflowerblue",
                       linestyle=":", alpha=1, label="1. Original Flux")

            for xvel in [cont_left[0], cont_left[1], cont_right[0],
                         cont_right[1]]:
                ax[0].axvline(xvel, color="darkblue", linestyle="-.", alpha=1)

            ax[0].plot(lsr_velocity[left[idx][0]:left[idx][1] + 1],
                       contadjspec.continuum[idx]
                       [left[idx][0]:left[idx][1] + 1],
                       color="darkblue", linestyle="-.", alpha=1,
                       label="2. Continuum Windows/Fits")
            ax[0].plot(lsr_velocity[right[idx][0]:right[idx][1] + 1],
                       contadjspec.continuum[idx]
                       [right[idx][0]:right[idx][1] + 1],
                       color="darkblue", linestyle="-.", alpha=1)

            flux_l = np.mean(raw_flux[left[idx][0]:left[idx][1]])
            flux_r = np.mean(raw_flux[right[idx][0]:right[idx][1]])
            mean_flux = np.mean([flux_l, flux_r])
            # norm = contadjspec.norm_flux[0] * 1e-14
            adj_delta = 0.25 * max_flux  # offset value, should make 0.25
            # scaling an option instead of hardcoded
            # filter out any flux values in the normalized flux that are less
            # than 0
            mask = np.where(contadjspec.norm_flux[idx] > 0)
            norm_filtered = contadjspec.norm_flux[idx][mask]
            norm_convolved = convolve(norm_filtered, kernel)
            lsr_filtered = lsr_velocity[mask]
            norm = norm_convolved * mean_flux + adj_delta
            # find the indices where the
            ax[0].step(lsr_filtered, norm, where="mid", data=None,
                       color="royalblue",
                       label=r"3. Normalized Flux ($\times {"r":.2e} + {"
                             r":.2e}$)".format(mean_flux, adj_delta))
            hflux_convolved = convolve(helper.flux[idx], kernel)
            hflux = hflux_convolved * mean_flux + adj_delta
            ax[0].step(helper.velocity[idx], hflux, where="mid", data=None,
                       color="darkblue", alpha=1,
                       label="4. Measurement Window")

            for xvel in [vel_min, vel_max]:
                ax[0].axvline(xvel, color="darkblue", linestyle="-", alpha=1)
                ax[1].axvline(xvel, color="darkblue", linestyle="-", alpha=1)

            acd = helper.acd[idx]["acd"]
            max_acd = max(acd)
            ax[1].set_ylim(0, 1.5 * max_acd)  # 1.5 is just a scaling factor to
            # make the plot prettier
            ax[1].bar(helper.velocity[idx], acd, width=helper.delta[idx],
                      color="maroon", label=r'5. Measured $N_{a}(v)$')

            ax[0].set_xlabel(r'\textbf{Velocity ($km/s$)}', fontsize=15)
            ax[0].xaxis.set_label_position('top')
            ax[0].xaxis.tick_top()
            ax[0].set_ylabel(r'\textbf{Flux ($ergs/cm^{2}/s/\AA$)}',
                             color="darkblue", fontsize=15)
            ax[1].set_xlabel(r'\textbf{Velocity ($km/s$)}', fontsize=15)
            ax[1].set_ylabel(r'\textbf{$N_{a}(v)$ ($ions/cm^{2}/(km/s)$)}',
                color="maroon", fontsize=15)

            ax[0].legend(loc="upper left", fontsize=12)
            ax[1].legend(loc="lower right", fontsize=12)

            ax[0].tick_params(axis="both", which="major", labelsize=12)
            ax[1].tick_params(axis="both", which="major", labelsize=12)

            ax[0].grid(alpha=0.4, linestyle="--", linewidth=1)
            ax[1].grid(alpha=0.4, linestyle="--", linewidth=1)

            if len(target) < 10:
                target_string = target
            else:
                target_string = target[0:10]

            outfile = OUTDIR + "{}_{}{}_visualization.png".format(
                                       target_string, ion, int(wv))

            if os.path.exists(outfile):
                os.remove(outfile)
                # it's getting corrupted if i don't do this?

            plt.savefig(os.path.join(outfile))
            LOGGER.info("Visualization {}/{} has been saved to {}.".format(
                idx + 1, len(helper.ions["ion"]), outfile))

            # once plt.show() is called it clears the figure so it has to be
            # last
            if show:
                plt.show()

        LOGGER.info("Visualization module complete. ")


# --------------------------------------------------------------------------- #
