#!/usr/bin/env python

""" Test how the predicted MSD changes as a function of the transition matrix using
the Markov State-Dependent model
"""

from LLC_Membranes.analysis.markov_state_dependent_dynamics import Chain
import numpy as np
import matplotlib.pyplot as plt

def transition_matrix(t):
	""" Make a 2D transition matrix with off diagonals 't' """

	T = np.zeros([2, 2])
	for i in range(2):
		for j in range(2):
			if i == j:
				T[i, j] = t
			else:
				T[i, j] = 1 - t
	return T

counts = 200000.  # total number of observations

state1 = [1.7, 0.0, 0.04]
state2 = [1.5, 0.0, 0.04]
transition = [1.4, 0.0, 0.04]

#H = np.array([0.1, 0.2, 0.4])[:, np.newaxis]
H = np.array([0.0, 0.0, 0.4])[:, np.newaxis]

emission_parameters = np.array([state1, state2, transition])

ts = [0.7, 0.85, 0.99]
#ts = [0.99, 0.7, 0.85]
#ts = [0.99]
for t in ts:

	count_matrix = counts * transition_matrix(t)
	chains = Chain(count_matrix, emission_parameters, hurst_parameters=H)	
	chains.generate_realizations(100, 1000, bound=1.0)
	chains.calculate_msd()
	chains.plot_msd(overlay=True, show=False, label='$T_{ii}$ =%.2f' % t)

plt.legend(fontsize=14)
plt.xlabel('Time (ns)', fontsize=14)
plt.ylabel('MSD (nm$^2$)', fontsize=14)
#plt.savefig('T_sensitivity.pdf')
plt.show()
