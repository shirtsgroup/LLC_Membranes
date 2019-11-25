#!/usr/bin/env python

""" Test how the predicted MSD changes as a function of the transition matrix using
the Markov State-Dependent model
"""

from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sol = ['URE', 'GCL', 'MET', 'ACH']
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
names2 = ['urea', 'ethylene\nglycol', 'methanol', 'acetic\nacid']
bar_width = 0.2
bar_locations = np.arange(1, 5)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

P = [0.57, 0.62, 0.62, 0.45]
PT = [0.40, 0.47, 0.44, 0.08]
PT_lambda = [0.0024, 0.0030, 0.0040, 0.0033]

params = np.array([P, PT, PT_lambda])

ax1.fill_between([1.05, 1.55], [0, 0], [0.7, 0.7], color='grey', alpha=0.3)
ax1.fill_between([2.05, 2.55], [0, 0], [0.7, 0.7], color='grey', alpha=0.3)
ax1.fill_between([3.05, 3.55], [0, 0], [0.7, 0.7], color='grey', alpha=0.3)
ax1.fill_between([4.05, 4.55], [0, 0], [0.7, 0.7], color='grey', alpha=0.3)

for i, res in enumerate(sol):

        ax1.bar(bar_locations[i] - bar_width, params[0, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        ax1.bar(bar_locations[i] - bar_width / 2 + 0.3, params[1, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        ax2.bar(bar_locations[i] - bar_width / 2 + 0.5, params[2, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7, hatch='//')

	#ax2.bar(bar_locations[i] + bar_width / 2, H[i], bar_width, hatch='//', edgecolor='black', color=colors[res], alpha=0.7)
	#plt.bar(bar_locations + i*bar_width - bar_width/2, heights, bar_width, label=names[res])

#ax1.text(0.65, 0.65, r'P($\alpha_d$)', fontsize=14)
#ax1.text(1.05, 0.65, r'P($\alpha_d, \lambda$)', fontsize=14)

hatch1 = mpatches.Patch(facecolor='white', label=r'$\alpha_d$', edgecolor='black')
hatch2 = mpatches.Patch(facecolor='white', hatch='//', label='$\lambda$', edgecolor='black')
plt.legend(handles=[hatch1, hatch2], fontsize=14, loc=0)
plt.xticks([1.05, 2.05, 3.05, 4.05], names2)
for xtick, res in zip(ax1.get_xticklabels(), sol):
    xtick.set_color(colors[res])

#ax1.set_xlabel('Solute', fontsize=14)
ax1.set_ylabel(r'$\alpha_d$', fontsize=14)
ax2.set_ylabel('$\lambda$', fontsize=14)
ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
ax1.set_ylim(0.0, 0.7)
#ax2.set_ylim(0.2, 0.45)
plt.tight_layout()
plt.savefig('1mode_AD_dwells.pdf')
plt.show()
