"""For debugging purposes only, instead of supplying the command line
arguments within the code."""

from spectrAOD.format_data import build_spectrum
from spectrAOD.measure_aod import parse, lsr_correct, continuum_fit, measure
from spectrAOD.visualization import Visualizer
from spectrAOD import SETTINGS

INPUTS = SETTINGS["inputs"]
DATADIR = INPUTS["datadir"]
PARAMETERS = SETTINGS["parameters"]
# set for debugging
ION = PARAMETERS["ion"]
ARGS = {
        "ion": ION,
        "instrument": PARAMETERS["instrument"],
        "filetype": PARAMETERS["filetype"],
        "vel_min": PARAMETERS["vel_min"],
        "vel_max": PARAMETERS["vel_min"],
        "grating": PARAMETERS["grating"],
        "redshift": PARAMETERS["redshift"]
        }


class Parameters:
    def __init__(self, ion, instrument, filetype, vel_min, vel_max, grating,
                 redshift):
        self.ion = ion
        self.instrument = instrument
        self.filetype = filetype
        self.vel_min = vel_min
        self.vel_max = vel_max
        self.grating = grating
        self.redshift = redshift


def main():
    # create the visualizer object
    visualizer = Visualizer()
    # argparse
    args, spectrum = parse()
    visualizer.set_target(spectrum.target)
    visualizer.set_raw_flux(spectrum.flux)
    # LSR correction
    spectrum = lsr_correct(args, spectrum)
    visualizer.set_raw_velocity(spectrum.raw_velocity)
    visualizer.set_lsr_velocity(spectrum.velocity)
    # continuum fit
    spectrum, left_indices, right_indices = continuum_fit(spectrum)
    visualizer.set_contadjspec(spectrum)
    visualizer.set_indices(left_indices, right_indices)
    # measure aod/acd/ew
    # set measurements back in spectrum object from helper object
    spectrum, helper = measure(args, spectrum)
    visualizer.set_helper(helper)
    # generate table
    spectrum.generate_table(args.vel_min, args.vel_max)
    visualizer.plot()
    return 0


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
