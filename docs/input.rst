.. _input-files:

Generate Gromacs Input Files
============================

*input.py* is a convenience script that will generate a GROMACS topology,
given an initial configuration, as well as .mdp files. Always check that
the files are accurate before using them.

*input.py* calls the classes defined on the following pages:

* :ref:`Generate Gromacs Topology File <gentop>`

* :ref:`Generate Gromacs .mdp Files <genmdp>`

.. _input-script:

==========================
Command Line Functionality
==========================

To use the script, run it from the command line with the following usage.

.. argparse::
   :filename: ../LLC_Membranes/setup/input.py
   :func: initialize
   :prog: input.py

