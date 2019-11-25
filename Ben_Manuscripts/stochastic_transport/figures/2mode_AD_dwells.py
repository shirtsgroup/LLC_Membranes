#!/usr/bin/env python

""" Test how the predicted MSD changes as a function of the transition matrix using
the Markov State-Dependent model
"""

from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#def align_yaxis(ax1, v1, ax2, v2):
#    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
#    _, y1 = ax1.transData.transform((0, v1))
#    _, y2 = ax2.transData.transform((0, v2))
#    inv = ax2.transData.inverted()
#    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
#    miny, maxy = ax2.get_ylim()
#    ax2.set_ylim(miny+dy, maxy+dy)

def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    adjust_yaxis(ax2,(y1-y2)/2,v2)
    adjust_yaxis(ax1,(y2-y1)/2,v1)

def adjust_yaxis(ax,ydif,v):
    """shift axis ax by ydiff, maintaining point v at the same location"""
    inv = ax.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, ydif))
    miny, maxy = ax.get_ylim()
    miny, maxy = miny - v, maxy - v
    if -miny>maxy or (-miny==maxy and dy > 0):
        nminy = miny
        nmaxy = miny*(maxy+dy)/(miny+dy)
    else:
        nmaxy = maxy
        nminy = maxy*(miny+dy)/(maxy+dy)
    ax.set_ylim(nminy+v, nmaxy+v)
sol = ['URE', 'GCL', 'MET', 'ACH']
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
names2 = ['urea', 'ethylene\nglycol', 'methanol', 'acetic\nacid']
bar_width = 0.2
bar_locations = np.arange(1, 5)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

# mode 1
P = [0.69, 0.69, 0.90, 0.58]
PT = [0.56, 0.62, 1.04, 0.41]
PT_lambda = [0.0037, 0.0026, 0.0006, 0.0026]

params1 = np.array([P, PT, PT_lambda])

# mode 2
P = [0.38, 0.48, 0.58, 0.33]
PT = [0.001, 0.06, 0.30, 0.001]
PT_lambda = [0.0027, 0.0049, 0.0054, 0.0021]

params2 = np.array([P, PT, PT_lambda])

ax1.fill_between([1.05, 1.55], [-1.1, -1.1], [1.1, 1.1], color='grey', alpha=0.3)
ax1.fill_between([2.05, 2.55], [-1.1, -1.1], [1.1, 1.1], color='grey', alpha=0.3)
ax1.fill_between([3.05, 3.55], [-1.1, -1.1], [1.1, 1.1], color='grey', alpha=0.3)
ax1.fill_between([4.05, 4.55], [-1.1, -1.1], [1.1, 1.1], color='grey', alpha=0.3)

for i, res in enumerate(sol):

        ax1.bar(bar_locations[i] - bar_width, params1[0, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        ax1.bar(bar_locations[i] - bar_width / 2 + 0.3, params1[1, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        ax2.bar(bar_locations[i] - bar_width / 2 + 0.5, params1[2, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7, hatch='//')

        ax1.bar(bar_locations[i] - bar_width, -params2[0, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        ax1.bar(bar_locations[i] - bar_width / 2 + 0.3, -params2[1, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        ax2.bar(bar_locations[i] - bar_width / 2 + 0.5, -params2[2, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7, hatch='//')

	#ax2.bar(bar_locations[i] + bar_width / 2, H[i], bar_width, hatch='//', edgecolor='black', color=colors[res], alpha=0.7)
	#plt.bar(bar_locations + i*bar_width - bar_width/2, heights, bar_width, label=names[res])

#ax1.text(0.65, 0.65, r'P($\alpha_d$)', fontsize=14)
#ax1.text(1.05, 0.65, r'P($\alpha_d, \lambda$)', fontsize=14)

hatch1 = mpatches.Patch(facecolor='white', label=r'$\alpha_d$', edgecolor='black')
hatch2 = mpatches.Patch(facecolor='white', hatch='//', label='$\lambda$', edgecolor='black')
plt.legend(handles=[hatch1, hatch2], fontsize=14, loc=0)
plt.xticks([1.05, 2.05, 3.05, 4.05], names2)
#ax1.set_yticklabels([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0], labels=[1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1.0])
ax1.set_yticklabels([1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1.0])
ax2.set_yticklabels([0.006, 0.004, 0.002, 0, 0.002, 0.004, 0.006])
for xtick, res in zip(ax1.get_xticklabels(), sol):
    xtick.set_color(colors[res])

#ax1.set_xlabel('Solute', fontsize=14)
ax1.set_ylabel(r'$\alpha_d$', fontsize=14)
ax2.set_ylabel('$\lambda$', fontsize=14)
ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
ax2.set_ylim(-0.006, 0.0045)
ax1.set_ylim(-0.65, 1.1)
align_yaxis(ax1, 0, ax2, 0)
plt.tight_layout()
plt.savefig('2mode_AD_dwells.pdf')
plt.show()
