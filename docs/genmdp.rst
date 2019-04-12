.. _genmdp:

===========================
Generate Gromacs .mdp Files
===========================

=======
Classes
=======

.. autoclass:: genmdp.SimulationMdp
   :members: __init__, write_em_mdp, write_npt_mdp, write_nvt_mdp, write_nve_mdp, add_pull_groups

==========================
Command Line Functionality
==========================

One can run this as a script from the command line as follows. However,
in most cases it is probably more convenient to use :ref:`input.py <input-files>`

.. argparse::
   :filename: ../LLC_Membranes/setup/genmdp.py
   :func: initialize
   :prog: genmdp.py

