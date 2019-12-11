#!/usr/bin/env python

""" Test how the predicted MSD changes as a function of the transition matrix using
the Markov State-Dependent model
"""

from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

root_dir = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
sol = ['URE', 'GCL', 'MET', 'ACH']
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
colors = {'sigma_normal': 'xkcd:blue', 'sigma_levy': 'xkcd:gold', 'alpha': 'xkcd:orangered'}
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
names2 = ['urea', 'ethylene\nglycol', 'methanol', 'acetic\nacid']
bar_width = 0.2
bar_locations = np.arange(1, 5)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

#sigma = [0.325459, 0.34, 0.34579, 0.2715]
#sigma2 = [0.21, 0.231, 0.22, 0.16]
#alpha = [1.84, 1.9144, 1.80, 1.72]

# did these by hand
#sigma_err = [0.004413, _, 0.00598, 0.004748]
#sigma2_err = [_, 0.00300, _, _]
#alpha_err = [_, 0.01771, _, _]

#params = np.array([sigma, sigma2, alpha])

# Old plotting format
#ax1.fill_between([1.05, 1.55], [0, 0], [0.37, 0.37], color='grey', alpha=0.3)
#ax1.fill_between([2.05, 2.55], [0, 0], [0.37, 0.37], color='grey', alpha=0.3)
#ax1.fill_between([3.05, 3.55], [0, 0], [0.37, 0.37], color='grey', alpha=0.3)
#ax1.fill_between([4.05, 4.55], [0, 0], [0.37, 0.37], color='grey', alpha=0.3)

for i, res in enumerate(sol):

        # old plotting format
        # ax1.bar(bar_locations[i] - bar_width, params[0, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        # ax1.bar(bar_locations[i] - bar_width / 2 + 0.3, params[1, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7)
        # ax2.bar(bar_locations[i] - bar_width / 2 + 0.5, params[2, i], bar_width, label=names[res], edgecolor='black', color=colors[res], alpha=0.7, hatch='//')

        # gaussian sigma and powerlaw alpha
        params = file_rw.load_object('%s/%s/10wt/forecast_%s_1statepowerlaw_gaussian.pl' % (root_dir, res, res))
        gauss_sigma = [] 
        for j in params.hop_parameters[0]:
            gauss_sigma.append(j[1])

        params = file_rw.load_object('%s/%s/10wt/forecast_%s_1state_powercut_levy.pl' % (root_dir, res, res))
        levy_sigma = [] 
        levy_alpha = []
        for j in params.hop_parameters[0]:
            levy_sigma.append(j[2])
            levy_alpha.append(j[0])

        ax1.bar(bar_locations[i] - bar_width, np.mean(gauss_sigma), bar_width, label=names[res], edgecolor='black', color=colors['sigma_normal'], alpha=1, yerr=np.std(gauss_sigma))
        ax1.bar(bar_locations[i], np.mean(levy_sigma), bar_width, label=names[res], edgecolor='black', color=colors['sigma_levy'], alpha=1, yerr=np.std(levy_sigma))
        ax2.bar(bar_locations[i] + bar_width, np.mean(levy_alpha), bar_width, label=names[res], edgecolor='black', color=colors['alpha'], alpha=1, yerr=np.std(levy_alpha))

	#ax2.bar(bar_locations[i] + bar_width / 2, H[i], bar_width, hatch='//', edgecolor='black', color=colors[res], alpha=0.7)
	#plt.bar(bar_locations + i*bar_width - bar_width/2, heights, bar_width, label=names[res])

#ax1.text(0.65, 0.65, r'P($\alpha_d$)', fontsize=14)
#ax1.text(1.05, 0.65, r'P($\alpha_d, \lambda$)', fontsize=14)

hatch1 = mpatches.Patch(facecolor=colors['sigma_normal'], label=r'$\mathcal{N}(\mathbf{\sigma})$', edgecolor='black')
hatch2 = mpatches.Patch(facecolor=colors['sigma_levy'], label=r'$L(\mathbf{\sigma}, \alpha_h)$', edgecolor='black')
hatch3 = mpatches.Patch(facecolor=colors['alpha'], label=r'$L(\sigma, \mathbf{\alpha_h})$', edgecolor='black')

plt.legend(handles=[hatch1, hatch2, hatch3], fontsize=14, loc='upper right')
plt.xticks([1.0, 2.0, 3.0, 4.0], names2)
#for xtick, res in zip(ax1.get_xticklabels(), sol):
#    xtick.set_color(colors[res])

ax1.set_ylabel('$\sigma$', fontsize=14)
ax2.set_ylabel(r'$\alpha_h$', fontsize=14)
ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
ax1.set_ylim(0.0, 0.40)
ax2.set_ylim(1.0, 2.15)
plt.tight_layout()
plt.savefig('1mode_AD_hops.pdf')
plt.show()
