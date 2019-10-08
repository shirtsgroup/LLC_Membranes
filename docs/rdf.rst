.. _rdf:

Radial Distribution Functions
=============================

For hexagonal phase LLC systems, we define the radial distribution 
function of a component as the number density of that component
as a function of the radial distance from the pore center.

The radial distance is calculated as the euclidean distance 
from the component coordinates to the pore center on the xy 
plane.

We normalize the counts of each component in each bin based on
the bin volume. In our case the bin volume is defined as:

.. math:: z\pi(r_2^2 - r_1^2)

where :math:`z` is the membrane thickness, and :math:`r_1` and :math:`r_2` are
the edges of adjacent bins where :math:`r_2 > r_1`.

.. figure:: images/radial_distribution_annulus.png
   :scale: 30 %
   :align: center

   Looking down onto the :math:`xy` plane of an atomistically rendered HII LLC membrane,
   a single bin used while calculating the radial distribution function is represented 
   by the area shaded between the two concentric black circles. In this case, the pore
   center is defined based on the average location of the blue spheres.

=======
Classes
=======

.. autoclass:: rdf.System
   :members: __init__, build_com, radial_distribution_function, build_spline, bootstrap, plot

=========
Functions
=========

.. automodule:: rdf 
   :members: grps

=========
Examples
=========

.. code-block:: python

   # imports
   from LLC_Membranes.analysis import rdf

.. code-block:: python

   # The most basic usage
   sys = rdf.System(gro, traj, residue, monomer)
   sys.radial_distribution_function()
   sys.bootstrap(200)  # 200 bootstrap trials
   sys.plot(show=True)

.. code-block:: python

   # Restrict residue to the center of mass of a group of constituent atoms
   atoms = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']  # name of carbon atoms in the monomer head group
   sys = rdf.System(gro, traj, residue, monomer)
   sys.radial_distribution_function()
   sys.bootstrap(200)  # 200 bootstrap trials
   sys.plot(show=True)

==========================
Command Line Functionality
==========================

.. argparse:: 
   :filename: ../LLC_Membranes/analysis/rdf.py
   :func: initialize
   :prog: rdf.py

=============
Example Usage
=============

.. code-block:: bash

   # Basic usage
   test.sh -t trajectory.trr -g coordinates.gro -r residue -b monomer

.. code-block:: bash

   # RDF of multiple residues or groups
   # The RDF of the center of mass of atoms C, C1, C2, C3, C4, C5 of residue_1 and the center
   # of mass of all atoms making up residue_2 will be plotted on top of each other
   test.sh -t trajectory.trr -g coordinates.gro -r residue_1 residue_2 -atoms C C1 C2 C3 C4 C5 -b monomer 


 
