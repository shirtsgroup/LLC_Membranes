.. _gentop:

==============================
Generate Gromacs Topology File
==============================

=======
Classes
=======

.. autoclass:: gentop.SystemTopology
   :members: __init__, write_top, add_residue, remove_residue

==========================
Command Line Functionality
==========================

One can run this as a script from the command line as follows. However, 
in most cases it is probably more convenient to use :ref:`input.py <input-files>`

.. argparse::
   :filename: ../LLC_Membranes/setup/gentop.py
   :func: initialize
   :prog: gentop.py
