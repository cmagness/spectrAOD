import os
import sys
import yaml

SETTINGS = None

if os.environ["SPECTRAOD_SETTINGS"]:
    with open(os.environ["SPECTRAOD_SETTINGS"]) as yamlfile:
        SETTINGS = yaml.safe_load(yamlfile)

else:
    files_in_current_dir = [f for f in os.listdir('.') if os.path.isfile(f)]
    for f in files_in_current_dir:
        if ".yaml" in f:
            with open(f) as yamlfile:
                SETTINGS = yaml.safe_load(yamlfile)
                # default SETTINGS file has empty values. if any are not
                # filled in, exit
                if not all(SETTINGS.values()):
                    sys.exit("All input and parameter fields in the "
                             "configuration file must be assigned. Please "
                             "refer to documentation on how to configure this "
                             "package. Exiting...")
                break

# trim whitespace in SETTINGS, make case insensitive

for parameter in SETTINGS["inputs"]:
    value = SETTINGS["inputs"][parameter]
    value = value.strip()

for parameter in SETTINGS["parameters"]:
    value = SETTINGS["parameters"][parameter]
    if type(value) == str:
        value = value.strip().upper()

if not SETTINGS["defaults"]["all_ions"] == "mini_ions.csv":
    value = SETTINGS["defaults"]["all_ions"]
    value = value.strip()

if not SETTINGS:
    sys.exit("No configuration file found. Please refer to documentation on "
             "how to configure this package. Exiting...")
