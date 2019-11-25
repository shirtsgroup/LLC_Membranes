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

sigma = [0.33, 0.34, 0.35, 0.27]
sigma2 = [0.21, 0.23, 0.22, 0.16]
alpha = [1.84, 1.92, 1.80, 1.72]

params = np.array([sigma, sigma2, alpha])

ax1.fill_between([1.05, 1.55], [0, 0], [0.37, 0.37], color='grey', alpha=0.3)
ax1.fill_between([2.05, 2.55], [0, 0], [0.37, 0.37], color='grey', alpha=0.3)
ax1.fill_between([3.05, 3.55], [0, 0], [0.37, 0.37], color='grey', alpha=0.3)
ax1.fill_between([4.05, 4.55], [0, 0], [0.37, 0.37], color='grey', alpha=0.3)

for i, res in enumerate(sol):

        ax1.bar(bar_locations[i] - bar_width, params[0, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        ax1.bar(bar_locations[i] - bar_width / 2 + 0.3, params[1, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        ax2.bar(bar_locations[i] - bar_width / 2 + 0.5, params[2, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7, hatch='//')

	#ax2.bar(bar_locations[i] + bar_width / 2, H[i], bar_width, hatch='//', edgecolor='black', color=colors[res], alpha=0.7)
	#plt.bar(bar_locations + i*bar_width - bar_width/2, heights, bar_width, label=names[res])

#ax1.text(0.65, 0.65, r'P($\alpha_d$)', fontsize=14)
#ax1.text(1.05, 0.65, r'P($\alpha_d, \lambda$)', fontsize=14)

hatch1 = mpatches.Patch(facecolor='white', label=r'$\alpha_h$', edgecolor='black', hatch='//')
hatch2 = mpatches.Patch(facecolor='white', label='$\sigma$', edgecolor='black')
plt.legend(handles=[hatch2, hatch1], fontsize=14, loc=0)
plt.xticks([1.05, 2.05, 3.05, 4.05], names2)
for xtick, res in zip(ax1.get_xticklabels(), sol):
    xtick.set_color(colors[res])

ax1.set_ylabel('$\sigma$', fontsize=14)
ax2.set_ylabel(r'$\alpha_h$', fontsize=14)
ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
ax1.set_ylim(0.0, 0.37)
ax2.set_ylim(1.0, 2.0)
plt.tight_layout()
plt.savefig('1mode_AD_hops.pdf')
plt.show()
