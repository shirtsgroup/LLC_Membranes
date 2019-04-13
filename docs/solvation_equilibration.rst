Solvate an Initial Configuration
================================

=======
Classes
=======

.. autoclass:: solvation_equilibration.System
   :members: __init__, query_database, equilibrate, calculate_pore_water, write_final_pore_configuration, place_water_tails, full_equilibration

========
Examples
========

.. code-block:: python

  # This pseudocode illustrates the commands that will build and solvate a system from scratch
  sys = System()
  
  while not sys.converged:
      sys.query_database()  # make a guess at the pore radius
      sys.equilibrate()  # Build system with guess radius, add water and make sure it is stable
      sys.calculate_pore_water() # Determine the total water in the pore. Stop loop if within tolerance

  sys.write_final_pore_configuration()  # write out coordinate file with correct no. of pore water
  sys.place_water_tails()  # put the appropriate amount of water in the tails
  sys.full_equilibration()  # Fully equilibrate the fully solvated system

==========================
Command Line Functionality
==========================

For convenience, a full equilibration can be achieved using a single python 
script which contains the above class and runs the sequence of commands in 
the example code block above. 

.. argparse::
   :filename: ../LLC_Membranes/setup/solvation_equilibration.py
   :func: initialize
   :prog: solvation_equilibration.py

