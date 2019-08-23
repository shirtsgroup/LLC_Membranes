.. _stats:

Statistical Analysis
====================

.. code-block:: python

    from LLC_Membranes.llclib import stats

This library contains a relatively sparse number of functions for calculating
some statistics. Most statistics are generated via bootstrapping for which the
implementation is specific to the problem, so they are not generalized here. 

=======
Classes
=======

.. autoclass:: stats.Cdf
   :members: __init__, cdf, random_sample, update_cdf

=========
Functions
=========

.. automodule:: stats
      :members: confidence_interval, outliers

