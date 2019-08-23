.. _timeseries:

Time Series Analysis
====================

.. code-block:: python

   from LLC_Membranes.llclib import timeseries

Commonly used functions for working with time series.

=======
Classes
=======

.. autoclass:: timeseries.VectorAutoRegression 
   :members: __init__

=========
Functions
=========

.. automodule:: timeseries
      :members: acf_slow, acf, autocov, msd_straightforward, msd, bootstrap_msd, step_autocorrelation, correlograms, switch_points, calculate_moving_average

