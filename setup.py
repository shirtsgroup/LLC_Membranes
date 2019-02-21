#! /usr/bin/env python

"""
LLC_Membranes repository
Used to analyze LLC membrane molecular dynamics trajectories
"""
from setuptools import setup

DOCLINES = __doc__.split("\n")

try:
    import mdtraj
except ModuleNotFoundError | ImportError:
    import subprocess
    command = 'conda install mdtraj'
    p = subprocess.Popen(command.split())
    p.wait()

setup(
    # Self-descriptive entries which should always be present
    name='LLC_Membranes',
    version = 0.1,
    description = 'Set up, simulate and analyze MD simulations of lyotropic liquid crystal membranes',
    author='Ben Coscia',
    license='BSD-3-Clause',
    install_requires=['numpy', 'mdtraj', 'matplotlib', 'ruptures', 'scipy', 'tqdm', 'fbm'], 
    packages = ['LLC_Membranes', 'LLC_Membranes.setup', 'LLC_Membranes.analysis', 'LLC_Membranes.llclib', 'LLC_Membranes.timeseries'],
)
