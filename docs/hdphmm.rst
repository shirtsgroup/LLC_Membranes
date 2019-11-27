.. _hdphmm:

The Infinite Hidden Markov Model
================================

Proper use of the python tools below requires some knowledge of the inner
workings of the infinite hidden Markov Model (iHMM). I recommend reading at
least the following background and ideally some of the listed references
before using the code.
 
:ref:`I already know how this works, skip to the real documentation. <ihmm-implementation>`

Hidden Markov Models
--------------------

A standard hidden Markov model (HMM) is used to infer the hidden (or latent)
states of a system given data observations (see Figure below). Methods exist
which allow one to determine the maximum likelihood latent state sequence from
the data and then map them back to a matrix which describes the probability of
transitioning between the latent states. The main drawback to this approach is
that the number of latent states must be known beforehand.

.. figure:: images/hmm.png
   :scale: 59 %
   :align: left

   (left) A process which can be described as an HMM produces observations (Y) that are emissions dependent on some sequence of hidden states (S). The type of emissions (illustrated as different shapes here) helps determine the likelihood of the hidden state which produced that emission. (center) There is an associated probability with transitioning between hidden states. In the above image, if a system started in state 1, it is possible to transition from state 1 to state 2, from state 1 to state 3 or or to remain in state 1. (right) The probabilities of transitions between states are typically stored in a probability transition matrix where each entry, P\ :sub:`ij` \, describes the probability of transitioning from state i to state j given that the system is already in state i.

There are a number of python tools that have already implemented HMMs. If you
already know the number of states, I recommend starting with them:

* `pomegranate <https://pomegranate.readthedocs.io/en/latest/>`_ : In addition to HMMs, implements a number of related probabilistic models.
* `PyEmma <http://emma-project.org/latest/>`_ : Tools for constructing Markov State Models using molecular dynamics data.
* `MSMBuilder <http://msmbuilder.org/3.8.0/>`_ : Similar to PyEmma

The infinite hidden Markov model
--------------------------------

The infinite hidden Markov model (iHMM), also called the hierarchical Dirichlet
process hidden markov model (HDPHMM), does not require the number of states to
be known. It uses a Dirichlet process prior in order to make guesses at the 
rows of the probability transition matrix with shape :math:`\infty \times \infty`.
Of course, in practice, we are not going to be doing math on an infinite sized 
matrix. Instead one can choose matrix dimensions that are big enough that all
hidden states will be identified. This may require some experimentation.

In practice, you will likely never encounter a situation where an unwieldy
number of states are found because sampling from a Dirichlet process benefits
from a "rich get richer"-type algorithm. By this, I mean that emissions tend to
get assigned to states that already exist. An easy way to visualize this  
  

.. _ihmm-implementation:

=======
Classes
=======

.. autoclass:: infinite_hidden_markov_model.InfiniteHMM
   :members: __init__, inference, summarize_results

======================
Command Line Interface
======================

.. argparse:: 
   :filename: ../LLC_Membranes/machine_learning/infinite_hidden_markov_model.py
   :func: initialize
   :prog: infinite_hidden_markov_model.py

