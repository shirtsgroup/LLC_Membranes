.. _random:

Generate Random Numbers
=======================

Methods for generating random numbers from a distribution based on random 
draws from a uniform distribution.

These were developed to debug a python script ported from MATLAB in which
many random numbers were generated. Given the same random seed, one can 
draw the same uniform random numbers in MATLAB and python (using numpy).
Unfortunately that doesn't work for other distribution due to differences
in implementation. 

=========
Functions
=========

.. automodule:: rand
      :members: randombeta, randombetavariate, randomexponential, randomgammaint, randomgamma, randomnormal, randomwishart, randombinomial, randomdirichlet 

