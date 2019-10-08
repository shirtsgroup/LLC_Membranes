.. _analysis:

Post-Simulation Trajectory Analysis
===================================

The scripts described below are used to perform all post-simulation analysis
and create the figures in our :ref:`papers <publications>`. Short descriptions
of each script are given below. Generally all of the scripts can be run from
the command line by passing arguments. Click on them to learn more about what
they do, the classes and functions that can be borrowed from them, and how to
run them from the command line.

Structural Analysis :code:`LLC_Membranes.analysis`
--------------------------------------------------

.. table::
   :widths: 30 70

   +----------------------------------------------+---------------------------------------------------------+
   |Script Name                                   |Description                                              |
   +==============================================+=========================================================+
   |:ref:`p2p.py <p2p>`                           |Calculate the average distance between pore centers      |
   +----------------------------------------------+---------------------------------------------------------+
   |:ref:`rdf.py <rdf>`                           |Calculate radial distribution function of a residue with |
   |                                              |respect to pore centers                                  |
   +----------------------------------------------+---------------------------------------------------------+
   |:ref:`correlation_function.py <correlation>`  |Calculate the correlation function of a residue or part  |
   |                                              |of a residue along an axis                               |
   +----------------------------------------------+---------------------------------------------------------+

Non-covalent Interactions :code:`LLC_Membranes.analysis`
--------------------------------------------------------

.. table::
   :widths: 30 70

   +-----------------------------------------------------+------------------------------------------------------------------------+
   |   Script Name                                       | Description                                                            |
   +=====================================================+========================================================================+
   |:ref:`hbonds.py <hbonds>`                            | Identify hydrogen bonds based on geometric distance and angle criteria |
   +-----------------------------------------------------+------------------------------------------------------------------------+
   |:ref:`coordination_number.py <coord-number>`         |                                                                        |
   +-----------------------------------------------------+------------------------------------------------------------------------+

Time Series Analysis :code:`LLC_Membranes.timeseries`
-----------------------------------------------------

.. table::
   :widths: 30 70

   +-----------------------------------------------------+------------------------------------------------------------------------+
   |   Script Name                                       | Description                                                            |
   +=====================================================+========================================================================+
   |:ref:`msd.py <msd>`                                  | Calculate the mean squared displacement and, optionally, the diffusion |
   |                                                     | coefficient of a residue.                                              |
   +-----------------------------------------------------+------------------------------------------------------------------------+
   |:ref:`coordinate_trace.py <z-trace>`                 | Track z-coordinate of a residue center of mass along with its radial   |
   |                                                     | position with respect to the pore center                               |
   +-----------------------------------------------------+------------------------------------------------------------------------+
   |:ref:`ctrwsim.py <ctrwsim>`                          | Simulate realizations of a continuous time random walk                 |
   +-----------------------------------------------------+------------------------------------------------------------------------+
   |:ref:`sfbm_parameters.py <sfbm-parameters>`          | Fit parameters assuming a solute undergoes subordinated fractional     |
   |                                                     | Brownian motion                                                        |
   +-----------------------------------------------------+------------------------------------------------------------------------+
   |:ref:`forecast_ctrw.py <forecast-ctrw>`              | Fit a correlated continuous time random walk to residue trajectory data|
   +-----------------------------------------------------+------------------------------------------------------------------------+
   |:ref:`fractional_levy_motion.py <flm>`               | Simulate realizations of fractional Levy motion                        |
   +-----------------------------------------------------+------------------------------------------------------------------------+
