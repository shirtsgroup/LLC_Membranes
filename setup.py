#! /usr/bin/env python

"""
LLC_Membranes repository
Used to analyze LLC membrane molecular dynamics trajectories
"""
from setuptools import setup
import versioneer

DOCLINES = __doc__.split("\n")

setup(
    # Self-descriptive entries which should always be present
    name='LLC_Membranes',
    author='Ben Coscia',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    install_requires=['numpy', 'mdtraj', 'matplotlib', 'ruptures', 'scipy', 'tqdm', 'fbm'], 

    # Which Python importable modules should be included when your package is installed
    # packages=['LLC_Membranes', 'LLC_Membranes.tests', 'LLC_Membranes.setup', 'LLC_Membranes.analysis', 'LLC_Membranes.llclib', 'LLC_Membranes.timeseries'],
    packages = ['tests', 'setup', 'analysis', 'llclib', 'timeseries'],
    # Optional include package data to ship with your package
    # Comment out this line to prevent the files from being packaged with your software
    # Extend/modify the list to include/exclude other items as need be
    # package_data={'fomms_integrate': ["data/*.dat"]
    #              },

    # Additional entries you may want simply uncomment the lines you want and fill in the data
    # author_email='me@place.org',      # Author email
    # url='http://www.my_package.com',  # Website
    # install_requires=[],              # Required packages, pulls from pip if needed; do not use for Conda deployment
    # platforms=['Linux',
    #            'Mac OS-X',
    #            'Unix',
    #            'Windows'],            # Valid platforms your code works on, adjust to your flavor
    # python_requires=">=3.5",          # Python version restrictions

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
