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

#sol = ['URE', 'GCL', 'MET', 'ACH']
sol = ['MET', 'ACH', 'URE', 'GCL']
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11'
bar_width = 0.18
bar_locations = np.arange(1, 10)
alpha = 0.7

# each of the following rows corresponds to a residue
alphah = [
    [1.79, 1.80, 1.88, 1.95, 1.34, 1.45, 1.61, 1.71, 1.42],
    [1.68, 1.75, 1.86, 1.91, 1.40, 1.52, 1.60, 1.74, 1.44],
    [1.56, 1.63, 1.80, 1.75, 1.28, 1.50, 1.20, 1.83, 1.45],
    [1.78, 1.88, 2.00, 2.00, 1.47, 1.70, 1.77, 2.00, 1.54],
    ]

for i, res in enumerate(sol):
	heights = alphah[i]
	plt.bar(bar_locations + (i - 1)*bar_width - bar_width/2, heights, bar_width, label=names[res], color=colors[res], edgecolor='black', alpha=alpha)

plt.legend(fontsize=14)
plt.xticks(ticks=bar_locations, labels=[1, 2, 3, 4, 5, 6, 7, 8, 'T'])
plt.xlabel('State', fontsize=14)
plt.ylabel(r'$\alpha_h$', fontsize=14)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.ylim(1, 2.25)
plt.yticks([1, 1.25, 1.5, 1.75, 2.0])
plt.savefig('alpha_v_state.pdf')
plt.show()
