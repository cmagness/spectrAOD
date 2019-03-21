import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

setuptools.setup(
    name="spectrAOD",
    version="0.0.1",
    author="Camellia Magness",
    author_email="cmagness@stsci.edu",
    description="This package is for measuring the apparent optical depth of spectra",
    long_description=long_description,
    long_description_content_type="text/markdown",
    # url="https://github.com/pypa/sampleproject",
    packages=setuptools.find_packages(),
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
)

