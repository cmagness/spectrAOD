SpectrAOD
------------
![spectraod_technique](spectraod_ddrf_final.png)

`spectrAOD` is a package for measuring the apparent optical depth, and thus,
 the apparent column density for spectral absorption features as outlined by
 Savage and Sembach in their 1991 paper. Many researchers have their own 
 code for performing these measurements but we set out to create an 
 open-source, well-maintained python package for people who didn't want to 
 reinvent the wheel. Currently, the package has a limited use case scenario 
 for Cosmic Origins Spectrograph data at a redshift of zero, but we hope to 
 expand the capabilities of the package in terms of what missions are 
 supported, the computation abilities, and visualization for these 
 calculations. Stay tuned!

Installation
------------

This package is registered on PyPI and available via `pip`. We also offer 
instructions on how to install from cloning the repository for the latest 
version, between published releases. You will need a working, and preferably current 
version of Anaconda.

##### Make a new environment
```
conda create --name <environment_name> python=3.5 <other packages>
```
Activate the new environment with:
```
conda activate <environment_name>
```
We recommend a short and simple name for the environment such as `spectraod`.

##### Install with `pip`

From the command line, in your new environment:
```
pip install spectrAOD==0.0.1
```

You can drop the version number and just use the name of the package if you 
would like the version most recently published.

##### Clone the repository and install it

This repository has a button near the top where you can click for the link 
to clone or download. Choose the https version unless you have set up an ssh
token for Github. 
Move into the directory that you would like this package to live in, then:
```
git clone https://github.com/cmagness/spectrAOD.git
cd spectrAOD
pip install .
```
Alternately, if you are having issues installing with `pip`, you can also 
use `python setup.py install`.

Using `spectrAOD`
-----------------

##### Configuring Settings

In this repository you will see a file called `sample_settings.yaml`. You 
will need a settings file to use this package that is of the same format. 
Copy the settings file and rename it as you please. For `spectrAOD` to find 
this file, you have two options:

###### Set an environment variable (Recommended)

We recommend creating an environment variable in your `.bashrc` or `
.bash_profile` to point to this file. To do so, open your `.bashrc` or `
.bash_profile` in a text editor (this will be a hidden file in your home 
directory if you are unfamiliar) and then add the line:
```
export SPECTRAOD_SETTINGS="/path/to/your/settings/file"
```
Save and close your `.bash_profile` and then activate these changes with:
```
source .bash_profile
```
Now `spectrAOD` will know where to look for your file. If you decide to move
it, just update the path.
 
###### Move your copy (Slightly faster)

If you don't want to mess with setting an environment variable, that is just
 fine. You can move your settings file to the directory you plan to run the 
 package from and `spectrAOD` will look for a `.yaml` file if no environment
  variable is set. Just be warned that it will look for _any_ `.yaml` file.
   
##### Running the package 

Once you've told `spectrAOD` where to find your settings file, be sure to 
actually populate it. Each parameter in the sample file has a comment with 
information about what should go into that file. Remove that comment and 
replace it with the default entries you want to use--don't worry, you can 
change the value in each run from the command line.

###### Setting a Target List

In the settings file you will notice one of the parameters asks for the path
 to your target list file. To perform the LSR correction `spectrAOD` needs a
  target list that has the RA and DEC. The target name needs to match the 
  name that is used in the file header if you are processing fits files. You
   can see the format for this file (it must also be a csv at this time) in 
   sample_targets.csv. Feel free to use that file to build your target list.

###### Command line arguments

To run the default measurements you've put in that file, from the command 
line, in any directory, enter:
```
measure <ion>
```
or, from within the package level:
```
python measure_aod.py <ion>
```

Where the ion is the name of the ion you wish to measure. Currently, you may
 use any of the ions in the `mini_ions.csv` list, or add your own in the 
 same format to that file.
 
To see the full list of parameter options, run:
```
measure --help
```

This will give you information about all the parameters you can change from 
the command line. By default, `spectrAOD` will use the values in your 
settings file but you can alter any of them from the default by adding the 
correct flag at the command line.

Say you perform a measurement of NV in the region from -100 to 100 km/s but 
then decide it may be advantageous to perform the measurement in the window 
from -100 to 150 km/s. You can do this without making changes to your 
settings file by running:
```
measure NV --vel_max 150
```

<!---

Contributing Code, Documentation, or Feedback
---------------------------------------------


3rd Party Libraries this package requires
-----------------------------------------


License
-------

---> 