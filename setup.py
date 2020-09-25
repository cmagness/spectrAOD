from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name="spectrAOD",
    version="1.1.0",
    author="Camellia Magness",
    author_email="cmagness@stsci.edu",
    description="This package is for measuring the apparent optical depth of "
                "spectra",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/cmagness/spectrAOD",
    keywords=['astronomy'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Topic :: Software Development :: Libraries :: Python Modules',
        ],
    python_requires='>=3.5',  # 3.5 and higher
    packages=find_packages(),
    install_requires=[
        'setuptools',
        'numpy',
        'astropy',
        'pandas',
        'argparse',
        'pyyaml'
        ],
    package_data={'spectrAOD': ['mini_ions.csv']},
    entry_points={'console_scripts': ['measure=spectrAOD.measure_aod:main',
                                      'batch=spectrAOD.batch_run:main']}
    # dependency_links=[],
    )
