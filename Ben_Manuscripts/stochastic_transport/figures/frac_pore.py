#!/usr/bin/env python

""" Test how the predicted MSD changes as a function of the transition matrix using
the Markov State-Dependent model
"""

from LLC_Membranes.analysis.markov_state_dependent_dynamics import States 
from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
names2 = ['methanol', 'acetic\nacid', 'urea', 'ethylene\nglycol']
H = [0.30, 0.34, 0.37, 0.40]
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11'
bar_width = 0.4
bar_locations = np.arange(1, 5)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

for i, res in enumerate(sol):
	states = file_rw.load_object('%s/%s/10wt/states.pl' %(path, res))
	w, v = np.linalg.eig(states.transition_matrix.T)
	heights = v[:, 0] / v[:, 0].sum()
	ax1.bar(bar_locations[i] - bar_width / 2, heights[:4].sum(), bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
	ax2.bar(bar_locations[i] + bar_width / 2, H[i], bar_width, hatch='//', edgecolor='black', color=colors[res], alpha=0.7)
	#plt.bar(bar_locations + i*bar_width - bar_width/2, heights, bar_width, label=names[res])

hatch1 = mpatches.Patch(facecolor='white', label='Fraction', edgecolor='black')
hatch2 = mpatches.Patch(facecolor='white', hatch='//', label='Hurst', edgecolor='black')
plt.legend(handles=[hatch1, hatch2], fontsize=14, loc=9)
plt.xticks(np.arange(1, 5), names2)
#ax1.set_xlabel('Solute', fontsize=14)
ax1.set_ylabel('Fraction of time spent in tails', fontsize=14)
ax2.set_ylabel('Hurst Parameter', fontsize=14)
ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
ax1.set_ylim(0.2, 1.0)
ax2.set_ylim(0.2, 0.45)
plt.tight_layout()
plt.savefig('frac_pore.pdf')
plt.show()
