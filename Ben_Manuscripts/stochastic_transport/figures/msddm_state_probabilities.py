#!/usr/bin/env python

""" Test how the predicted MSD changes as a function of the transition matrix using
the Markov State-Dependent model
"""

from LLC_Membranes.analysis.markov_state_dependent_dynamics import States 
from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt
import tqdm

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
nboot = 200 

hatch1 = '///'
hatch2 = '...'
hatches = [hatch1, hatch1, hatch1, hatch1, hatch2, hatch2, hatch2, hatch2, None]

#for i, res in enumerate(sol):
#    states = file_rw.load_object('%s/%s/10wt/states.pl' %(path, res))
#    w, v = np.linalg.eig(states.transition_matrix.T)
#    heights = v[:, 0] / v[:, 0].sum()
#    for j in range(9):
#        plt.bar(bar_locations[j] + (i - 1)*bar_width - bar_width/2, heights[j], bar_width, label=names[res], color=colors[res], edgecolor='black', alpha=opacity, hatch=hatches[j])

for i, res in enumerate(sol):

        states = file_rw.load_object('%s/%s/10wt/states.pl' %(path, res))

        heights = np.zeros([nboot, 8])
        for b in tqdm.tqdm(range(nboot)):

            traj = np.random.choice(24, size=24, replace=True)
            tmatrix = np.zeros([8, 8])
            count_matrix = np.zeros([8, 8])
            state_sequence = states.state_sequence[:, traj]

            for t in range(1, states.nT):
                transitioned_from = state_sequence[t - 1, :]
                transitioned_to = state_sequence[t, :]
                for pair in zip(transitioned_from, transitioned_to):
                    count_matrix[pair[0], pair[1]] += 1 
             
            tmatrix = (count_matrix.T / count_matrix.sum(axis=1)).T
            w, v = np.linalg.eig(tmatrix.T)
            heights[b, :] = v[:, 0] / v[:, 0].sum()

#	w, v = np.linalg.eig(states.transition_matrix.T)
#	heights = v[:, 0] / v[:, 0].sum()
	#print(sum(heights[:4]))
	#print(w[0])  # make sure this is 1
	#print(v[:, 0] / v[:, 0].sum())
        for j in range(8):
            plt.bar(bar_locations[j] + (i - 1)*bar_width - bar_width/2, heights[:, j].mean(), bar_width, color=colors[res], edgecolor='black', alpha=opacity, hatch=hatches[j], yerr=heights[:, j].std())

        #plt.bar(bar_locations + (i - 1)*bar_width - bar_width/2, heights.mean(axis=0), bar_width, label=names[res], color=colors[res], edgecolor='black', alpha=opacity, yerr=heights.std(axis=0))


import matplotlib.patches as mpatches
hatch1 = mpatches.Patch(facecolor='white', label='In Tails', edgecolor='black', hatch=hatch1)
hatch2 = mpatches.Patch(facecolor='white', label='In Pores', edgecolor='black', hatch=hatch2)
patch1 = mpatches.Patch(facecolor=colors['URE'], label=names['URE'], edgecolor='black', alpha=opacity)
patch2 = mpatches.Patch(facecolor=colors['GCL'], label=names['GCL'], edgecolor='black', alpha=opacity)
patch3 = mpatches.Patch(facecolor=colors['MET'], label=names['MET'], edgecolor='black', alpha=opacity)
patch4 = mpatches.Patch(facecolor=colors['ACH'], label=names['ACH'], edgecolor='black', alpha=opacity)
patch5 = mpatches.Patch(facecolor='w', label='', edgecolor='w', alpha=opacity)
patch6 = mpatches.Patch(facecolor='w', label='', edgecolor='w', alpha=opacity)

#plt.rc('text', usetex=True)
labels = ['1\n$^{t}$', '2\n$^{(t/h)}$', '3\n$^{(t/a)}$', '4\n$^{(t/h/a)}$', '5\n$^{(p)}$', '6\n$^{(p/h)}$', '7\n$^{(p/a)}$', '8\n$^{(p/h/a)}$', 'T']

plt.legend(fontsize=14, handles=[patch1, patch2, patch3, patch4, hatch1, hatch2, patch5, patch6], ncol=2, columnspacing=1)
plt.xticks(ticks=bar_locations, labels=labels)
plt.xlabel('State', fontsize=14)
plt.ylabel('Probability of Occupation', fontsize=14)
plt.tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('state_probabilities.pdf')
plt.show()
