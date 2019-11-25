#!/usr/bin/env python

""" Test how the predicted MSD changes as a function of the transition matrix using
the Markov State-Dependent model
"""

from LLC_Membranes.analysis.markov_state_dependent_dynamics import States 
from LLC_Membranes.llclib import file_rw
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

sol = ['URE', 'GCL', 'MET', 'ACH']
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11'
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
bar_width = 0.18
bar_locations = np.arange(1, 9)
opacity = 0.7

for i, res in enumerate(sol):
	states = file_rw.load_object('%s/%s/10wt/states.pl' %(path, res))
	w, v = np.linalg.eig(states.transition_matrix.T)
	heights = v[:, 0] / v[:, 0].sum()
	#print(sum(heights[:4]))
	#print(w[0])  # make sure this is 1
	#print(v[:, 0] / v[:, 0].sum())
	plt.bar(bar_locations + (i - 1)*bar_width - bar_width/2, heights, bar_width, label=names[res], color=colors[res], edgecolor='black', alpha=opacity)

plt.legend(fontsize=14)
plt.xlabel('State', fontsize=14)
plt.ylabel('Probability of Occupation', fontsize=14)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('state_probabilities.pdf')
plt.show()
