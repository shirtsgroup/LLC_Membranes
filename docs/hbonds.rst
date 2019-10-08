.. _hbonds:

Identify Hydrogen Bonds
=======================

Use the functionality of hbonds.py in order to identify hydrogen bonds in a
molecular dynamics trajectory.

Attempts to describe a hydrogen bond in the context of molecular simulations
has yielded a number of definitions with no true consensus
\cite{prada-gracia_quest_2013} especially since the geometry of hydrogen bonds
has some dependence on the system being studied. Luzar and Chandler proposed
the geometric criterion such that a hydrogen bond exists if the distance
between the donor, D, and acceptor, A, atoms is less than 3.5 \AA~and the angle
formed by D--H$\cdot\cdot\cdot$A is less than 30\degree.
\cite{luzar_effect_1996} The definition of Luzar and Chandler is easily
visualized for trajectories using the \texttt{hbonds} representation of the
Visual Molecular Dynamics (VMD) software package which allows us to directly
check the validity of identified hydrogen bonds.


Some definitions:

* donor (D): atom covalently bonded to hydrogen
* acceptor (A) : atom which 'accepts' hydrogen bond

A diagram of an hbond:

    D--H - - A

Criterion:

* Distance between D and A below some distance
* Angle between DHA less than some cut-off

Note: You will need to properly annotate any residue that you think might
participate in a hydrogen bond. See :ref:`annotation-table`.

=======
Classes
=======

.. autoclass:: hbonds.System
   :members: __init__, set_eligible, identify_hbonds, plot_hbonds

=========
Examples
=========

.. code-block:: python

   # imports
   from LLC_Membranes.analysis import hbonds 

.. code-block:: python

   # Initialize system
   sys = hbonds.System(traj, gro)
   
   # Define residues involved in hbonds of interest
   residues = ['HII', 'HOH'] # hydrogen bonds between LLC monomer and water
   
   # Restrict system to hbonding atoms that are a part of HII and HOH 
   for i, r, in enumerate(residues):
       sys.set_eligible(r, 'all')  # restrict to all hbonding atoms of each residue

   # Apply geometric criteria to identify hydrogen bonds
   distance = 0.35  # distance cut-off 
   angle = 30
   sys.identify_hbonds(distance, angle)

   # Plot total hydrogen bonds as a function of time
   sys.plot_hbonds()

==========================
Command Line Functionality
==========================

.. argparse:: 
   :filename: ../LLC_Membranes/analysis/hbonds.py
   :func: initialize
   :prog: hbonds.py

=============
Example Usage
=============

.. code-block:: bash

  # hydrogen bonds between atoms O3 and O4 of the monomer HII with all atoms of water
  hbonds.py -t trajectory.xtc -g coordinates.gro -r HII HOH -a O3 O4 -a all

