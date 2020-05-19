"""For debugging purposes only, instead of supplying the command line
arguments within the code."""

from spectrAOD.format_data import build_spectrum
from spectrAOD.measure_aod import lsr_correct, continuum_fit, measure
from spectrAOD import SETTINGS

INPUTS = SETTINGS["inputs"]
DATADIR = INPUTS["datadir"]
PARAMETERS = SETTINGS["parameters"]
# set for debugging
ION = "NiII"
ARGS = {
        "ion": ION,
        "instrument": PARAMETERS["instrument"],
        "filetype": PARAMETERS["filetype"],
        "vel_min": PARAMETERS["vel_min"],
        "vel_max": PARAMETERS["vel_min"],
        "grating": PARAMETERS["grating"]
        }


class Parameters:
    def __init__(self, ion, instrument, filetype, vel_min, vel_max, grating):
        self.ion = ion
        self.instrument = instrument
        self.filetype = filetype
        self.vel_min = vel_min
        self.vel_max = vel_max
        self.grating = grating


def main():
    # argparse
    args = Parameters(**ARGS)
    spectrum = build_spectrum(DATADIR, args.instrument.upper(),
                              args.filetype.upper(), args.grating.upper())
    # LSR correction
    spectrum = lsr_correct(args, spectrum)
    # continuum fit
    spectrum = continuum_fit(spectrum)
    # measure aod/acd/ew
    # set measurements back in spectrum object from helper object
    spectrum = measure(args, spectrum)
    # generate table
    spectrum.generate_table(args.vel_min, args.vel_max)


# -----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
