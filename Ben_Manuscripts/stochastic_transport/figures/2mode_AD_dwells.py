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

root_dir = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
sol = ['URE', 'GCL', 'MET', 'ACH']
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
colors = {'powerlaw_alpha': 'xkcd:blue', 'powercut_alpha': 'xkcd:gold', 'powercut_lambda': 'xkcd:orangered'}
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
names2 = ['urea', 'ethylene\nglycol', 'methanol', 'acetic\nacid']
bar_width = 0.2
bar_locations = np.arange(1, 5)
fontsize=16
opacity = 1

fig, ax1 = plt.subplots(figsize=(8,6))
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

#ax1.fill_between([1.05, 1.55], [-1.1, -1.1], [1.1, 1.1], color='grey', alpha=0.3)
#ax1.fill_between([2.05, 2.55], [-1.1, -1.1], [1.1, 1.1], color='grey', alpha=0.3)
#ax1.fill_between([3.05, 3.55], [-1.1, -1.1], [1.1, 1.1], color='grey', alpha=0.3)
#ax1.fill_between([4.05, 4.55], [-1.1, -1.1], [1.1, 1.1], color='grey', alpha=0.3)

for i, res in enumerate(sol):

        params = file_rw.load_object('%s/%s/10wt/forecast_%s_2state_powerlaw_gaussian.pl' % (root_dir, res, res))
        powerlaw_alpha1 = []
        for j in params.dwell_parameters[0]:
            powerlaw_alpha1.append(j)
        powerlaw_alpha2 = []
        for j in params.dwell_parameters[1]:
            powerlaw_alpha2.append(j)

        params = file_rw.load_object('%s/%s/10wt/forecast_%s_2state_powercut_levy.pl' % (root_dir, res, res))
        powercut_alpha1 = []
        powercut_lambda1 = []
        for j in params.dwell_parameters[0]:
            powercut_alpha1.append(j[0])
            powercut_lambda1.append(j[1])
        powercut_alpha2 = []
        powercut_lambda2 = []
        for j in params.dwell_parameters[1]:
            powercut_alpha2.append(j[0])
            powercut_lambda2.append(j[1])

        ax1.bar(bar_locations[i] - bar_width, -np.mean(powerlaw_alpha1), bar_width, label=names[res], edgecolor='black', color=colors['powerlaw_alpha'], alpha=opacity, yerr=np.std(powerlaw_alpha1))
        ax1.bar(bar_locations[i], -np.mean(powercut_alpha1), bar_width, label=names[res], edgecolor='black', color=colors['powercut_alpha'], alpha=opacity, yerr=np.std(powercut_alpha1))
        ax2.bar(bar_locations[i] + bar_width, -np.mean(powercut_lambda1), bar_width, label=names[res], edgecolor='black', color=colors['powercut_lambda'], alpha=opacity, yerr=np.std(powercut_lambda1))

        ax1.bar(bar_locations[i] - bar_width, np.mean(powerlaw_alpha2), bar_width, label=names[res], edgecolor='black', color=colors['powerlaw_alpha'], alpha=opacity, yerr=np.std(powerlaw_alpha2))
        ax1.bar(bar_locations[i], np.mean(powercut_alpha2), bar_width, label=names[res], edgecolor='black', color=colors['powercut_alpha'], alpha=opacity, yerr=np.std(powercut_alpha2))
        ax2.bar(bar_locations[i] + bar_width, np.mean(powercut_lambda2), bar_width, label=names[res], edgecolor='black', color=colors['powercut_lambda'], alpha=opacity, yerr=np.std(powercut_lambda2))

        #ax1.bar(bar_locations[i] - bar_width, params1[0, i], bar_width, label=names[res], edgecolor='black', color=colors['powerlaw_alpha'], alpha=opacity)
        #ax1.bar(bar_locations[i], params1[1, i], bar_width, label=names[res], edgecolor='black', color=colors['powercut_alpha'], alpha=opacity)
        #ax2.bar(bar_locations[i] + bar_width, params1[2, i], bar_width, label=names[res], edgecolor='black', color=colors['powercut_lambda'], alpha=opacity)

        #ax1.bar(bar_locations[i] - bar_width, -params2[0, i], bar_width, label=names[res], edgecolor='black', color=colors['powerlaw_alpha'], alpha=opacity)
        #ax1.bar(bar_locations[i], -params2[1, i], bar_width, label=names[res], edgecolor='black', color=colors['powercut_alpha'], alpha=opacity)
        #ax2.bar(bar_locations[i] + bar_width, -params2[2, i], bar_width, label=names[res], edgecolor='black', color=colors['powercut_lambda'], alpha=opacity)

        #ax1.bar(bar_locations[i] - bar_width, params1[0, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        #ax1.bar(bar_locations[i] - bar_width / 2 + 0.3, params1[1, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        #ax2.bar(bar_locations[i] - bar_width / 2 + 0.5, params1[2, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7, hatch='//')

        #ax1.bar(bar_locations[i] - bar_width, -params2[0, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        #ax1.bar(bar_locations[i] - bar_width / 2 + 0.3, -params2[1, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        #ax2.bar(bar_locations[i] - bar_width / 2 + 0.5, -params2[2, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7, hatch='//')

	#ax2.bar(bar_locations[i] + bar_width / 2, H[i], bar_width, hatch='//', edgecolor='black', color=colors[res], alpha=0.7)
	#plt.bar(bar_locations + i*bar_width - bar_width/2, heights, bar_width, label=names[res])

#ax1.text(0.65, 0.65, r'P($\alpha_d$)', fontsize=14)
#ax1.text(1.05, 0.65, r'P($\alpha_d, \lambda$)', fontsize=14)

hatch1 = mpatches.Patch(facecolor=colors['powerlaw_alpha'], label=r'$P(\mathbf{\alpha_d})$', edgecolor='black')
hatch2 = mpatches.Patch(facecolor=colors['powercut_alpha'], label=r'$P_T(\mathbf{\alpha_d}, \lambda)$', edgecolor='black')
hatch3 = mpatches.Patch(facecolor=colors['powercut_lambda'], label=r'$P_T(\alpha_d, \mathbf{\lambda}$', edgecolor='black')

plt.legend(handles=[hatch1, hatch2, hatch3], fontsize=fontsize, loc=0)
plt.xticks([1.0, 2.0, 3.0, 4.0], names2)
#ax1.set_yticklabels([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0], labels=[1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1.0])
ax1.set_yticklabels([1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1.0])
ax2.set_yticklabels([0.006, 0.004, 0.002, 0, 0.002, 0.004, 0.006])
#for xtick, res in zip(ax1.get_xticklabels(), sol):
#    xtick.set_color(colors[res])

ax1.text(2.6, 1, 'Pores', verticalalignment='center', horizontalalignment='center', fontsize=18, fontweight='bold')#, bbox=props, fontweight='bold')
ax1.text(2.6, -0.9, 'Tails', verticalalignment='center', horizontalalignment='center', fontsize=18, fontweight='bold')#, bbox=props, fontweight='bold')

#ax1.set_xlabel('Solute', fontsize=14)
ax1.set_ylabel(r'$\alpha_d$', fontsize=fontsize)
ax2.set_ylabel('$\lambda$', fontsize=fontsize)
ax1.tick_params(labelsize=fontsize)
ax2.tick_params(labelsize=fontsize)
ax1.plot([0.6, 4.6], [0, 0], lw=2, color='black')
ax1.set_xlim(0.6, 4.6)
ax2.set_ylim(-0.006, 0.0045)
ax1.set_ylim(-0.65, 1.1)
align_yaxis(ax1, 0, ax2, 0)
plt.tight_layout()
plt.savefig('2mode_AD_dwells.pdf')
plt.show()
