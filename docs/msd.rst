.. _msd:

Calculate Mean Squared Displacement
===================================

A particle's mean squared displacement (MSD) tells you the squared distance it
moves over time. By analyzing the shape of the MSD curve, we can gain
mechanistic insight. The general form of the MSD curve is:

.. math::

   \langle x^2(t) \rangle = K_{\alpha}t^{\alpha}

where :math:`K_{\alpha}` is the generalized diffusion coefficient. In cases
were :math:`\alpha=1`, the MSD is linear and :math:`K_{\alpha}` is equal to
the standard diffusion coefficient measured for particles undergoing
Brownian motion. Albert Einstein was the first to prove that the MSD of a
Brownian particle is linear and that its slope is proportional to the diffusion
coefficient :math:`D`:

.. math::

   \frac{d MSD}{dt} \propto 2NDt

where :math:`N` is the number of spatial dimensions over which the MSD
was measured. When :math:`\alpha > 1` diffusion occurs faster than expected
and is termed superdiffusion. When :math:`\alpha < 1`, diffusion is slower
than expected and called subdiffusion.

There are two ways one can calculate MSD. The first and simplest way is 
the ensemble MSD:

.. math:: 

   \langle x^2(t) \rangle = \langle x(t) - x(0) \rangle^2

The ensemble MSD measures the distance traveled from a particle's initial
position. The second way measures and averages the particle's displacement over
all possible time lags, :math:`\tau`, and is called the time-averaged MSD:

.. math::

   \overline{x^2(\tau)} = \dfrac{1}{T - \tau}\int_{0}^{T - \tau} (x(t + \tau) - z(t))^2 dt

where :math:`T` is the length of the trajectory. The time-averaged MSD is
more statistically robust than the ensemble average. For ergodic systems,
both methods of calculation will yield the same result. However, with 
less data, the time-averaged MSD will give tighter error bars.

:code:`msd.py` can calculate both types of MSDs for a particle. If you believe
the particle undergoes Brownian motion, then you can also measure the diffusion
constant. The script can be run from the :ref:`command line <msd-command-line>`
or its :ref:`classes <msd-classes>` can be imported for use in other analysis
scripts. Both usages are explained below.

.. _msd-command-line:

======================
Command Line Interface
======================

.. argparse::
   :filename: ../LLC_Membranes/timeseries/msd.py
   :func: initialize
   :prog: msd.py

.. _msd-classes:

========
Classes
========

.. autoclass:: msd.Diffusivity
   :members: __init__, restrict_to_pore, calculate, step_autocorrelation, fit_linear, fit_power_law, stationarity, bootstrap_power_law, bootstrap, plot, plot_power_law, plot_autocorrelation
