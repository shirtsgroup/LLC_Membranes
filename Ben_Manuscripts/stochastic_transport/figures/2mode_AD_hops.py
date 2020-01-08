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
colors = {'URE':'xkcd:magenta', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
colors = {'sigma_normal': 'xkcd:magenta', 'sigma_levy': 'xkcd:orange', 'alpha': 'xkcd:blue'}
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
names2 = ['urea', 'ethylene\nglycol', 'methanol', 'acetic\nacid']
bar_width = 0.2
bar_locations = np.arange(1, 5)
fontsize=16
opacity=1

fig, ax1 = plt.subplots(figsize=(8, 6))
ax2 = ax1.twinx()

# mode 1
sigma = [0.35, 0.38, 0.45, 0.32]
sigmaL = [0.24, 0.26, 0.31, 0.21]
alpha = [1.91, 1.99, 1.97, 1.91]

params1 = np.array([sigma, sigmaL, alpha])

# mode 2

sigma = [0.24, 0.23, 0.32, 0.17]
sigmaL = [0.12, 0.15, 0.20, 0.09]
alpha = [1.50, 1.90, 1.85, 1.50]

params2 = np.array([sigma, sigmaL, alpha])

lowerlim = [-1.5, -1.5]
upperlim = [1.5, 1.5]

#ax1.fill_between([1.05, 1.55], lowerlim, upperlim, color='grey', alpha=0.3)
#ax1.fill_between([2.05, 2.55], lowerlim, upperlim, color='grey', alpha=0.3)
#ax1.fill_between([3.05, 3.55], lowerlim, upperlim, color='grey', alpha=0.3)
#ax1.fill_between([4.05, 4.55], lowerlim, upperlim, color='grey', alpha=0.3)

for i, res in enumerate(sol):

        params = file_rw.load_object('%s/%s/10wt/forecast_%s_2state_powerlaw_gaussian.pl' % (root_dir, res, res))
        gauss_sigma1 = []
        gauss_sigma2 = []
        for j in params.hop_parameters[0]:
            gauss_sigma1.append(j[1])
        for j in params.hop_parameters[1]:
            gauss_sigma2.append(j[1])

        params = file_rw.load_object('%s/%s/10wt/forecast_%s_2state_powercut_levy.pl' % (root_dir, res, res))
        levy_sigma1 = []
        levy_alpha1 = []
        levy_sigma2 = []
        levy_alpha2 = []
        for j in params.hop_parameters[0]:
            levy_sigma1.append(j[2])
            levy_alpha1.append(j[0])
        for j in params.hop_parameters[1]:
            levy_sigma2.append(j[2])
            levy_alpha2.append(j[0])

        ax1.bar(bar_locations[i] - bar_width, -np.mean(gauss_sigma1), bar_width, label=names[res], edgecolor='black', color=colors['sigma_normal'], alpha=opacity, yerr=np.std(gauss_sigma1))
        ax1.bar(bar_locations[i], -np.mean(levy_sigma1), bar_width, label=names[res], edgecolor='black', color=colors['sigma_levy'], alpha=opacity, yerr=np.std(levy_sigma1))
        ax2.bar(bar_locations[i] + bar_width, -np.mean(levy_alpha1), bar_width, label=names[res], edgecolor='black', color=colors['alpha'], alpha=opacity, yerr=np.std(levy_alpha1))

        ax1.bar(bar_locations[i] - bar_width, np.mean(gauss_sigma2), bar_width, label=names[res], edgecolor='black', color=colors['sigma_normal'], alpha=opacity, yerr=np.std(gauss_sigma2))
        ax1.bar(bar_locations[i], np.mean(levy_sigma2), bar_width, label=names[res], edgecolor='black', color=colors['sigma_levy'], alpha=opacity, yerr=np.std(levy_sigma2))
        ax2.bar(bar_locations[i] + bar_width, np.mean(levy_alpha2), bar_width, label=names[res], edgecolor='black', color=colors['alpha'], alpha=opacity, yerr=np.std(levy_alpha2))

        #ax1.bar(bar_locations[i] - bar_width, np.mean(gauss_sigma1), bar_width, label=names[res], edgecolor='black', color=colors['sigma_normal'], alpha=opacity, yerr=np.std(gauss_sigma1))
        #ax1.bar(bar_locations[i], np.mean(levy_sigma1), bar_width, label=names[res], edgecolor='black', color=colors['sigma_levy'], alpha=opacity)
        #ax2.bar(bar_locations[i] + bar_width, params1[2, i], bar_width, label=names[res], edgecolor='black', color=colors['alpha'], alpha=opacity)

        #ax1.bar(bar_locations[i] - bar_width, -params2[0, i], bar_width, label=names[res], edgecolor='black', color=colors['sigma_normal'], alpha=opacity)
        #ax1.bar(bar_locations[i], -params2[1, i], bar_width, label=names[res], edgecolor='black', color=colors['sigma_levy'], alpha=opacity)
        #ax2.bar(bar_locations[i] + bar_width, -params2[2, i], bar_width, label=names[res], edgecolor='black', color=colors['alpha'], alpha=opacity)

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

hatch1 = mpatches.Patch(facecolor=colors['sigma_normal'], label=r'$\mathcal{N}(\mathbf{\sigma})$', edgecolor='black')
hatch2 = mpatches.Patch(facecolor=colors['sigma_levy'], label=r'$L(\mathbf{\sigma}, \alpha_h)$', edgecolor='black')
hatch3 = mpatches.Patch(facecolor=colors['alpha'], label=r'$L(\sigma, \mathbf{\alpha_h})$', edgecolor='black')

plt.legend(handles=[hatch1, hatch2, hatch3], fontsize=14, loc='lower left')
plt.xticks([1.0, 2.0, 3.0, 4.0], names2)
#ax1.set_yticklabels([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0], labels=[1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1.0])
#ax1.set_yticklabels([1, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1.0])
#ax2.set_yticklabels([0.006, 0.004, 0.002, 0, 0.002, 0.004, 0.006])
#for xtick, res in zip(ax1.get_xticklabels(), sol):
#    xtick.set_color(colors[res])

props = dict(boxstyle='square', facecolor='grey', alpha=0.3, lw=0, fill=False)
ax1.text(2.6, 1.2, 'Pores', verticalalignment='center', horizontalalignment='center', fontsize=18, fontweight='bold')#, bbox=props, fontweight='bold')
ax1.text(2.6, -1.2, 'Tails', verticalalignment='center', horizontalalignment='center', fontsize=18, fontweight='bold')#, bbox=props, fontweight='bold')

ax1.set_yticklabels([1.5, 1, 0.5, 0, 0.5, 1])#, 0.75, 0.5, 0.25, 0, 0.25, 0.5, 0.75, 1.0])
ax2.set_yticklabels([3, 2, 1, 0, 1, 2])#0.006, 0.004, 0.002, 0, 0.002, 0.004, 0.006])

ax1.plot([0, 5], [0, 0], color='black', lw=2)
#ax1.set_xlabel('Solute', fontsize=14)
ax1.set_ylabel('$\sigma$', fontsize=fontsize)
ax2.set_ylabel(r'$\alpha_h$', fontsize=fontsize)
ax1.tick_params(labelsize=fontsize)
ax2.tick_params(labelsize=fontsize)
ax2.set_ylim(-2.5, 2.5)
ax1.set_ylim(-1.3, 1.3)
ax1.set_xlim(0.6, 4.6)
align_yaxis(ax1, 0, ax2, 0)
plt.tight_layout()
plt.savefig('2mode_AD_hops.pdf')
plt.show()
