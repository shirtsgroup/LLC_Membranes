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
colors = {'powerlaw_alpha': 'xkcd:blue', 'powercut_alpha': 'xkcd:gold', 'powercut_lambda': 'xkcd:orangered'}
names = {'URE': 'urea', 'GCL': 'ethylene glycol', 'MET': 'methanol', 'ACH': 'acetic acid'}
names2 = ['urea', 'ethylene\nglycol', 'methanol', 'acetic\nacid']
bar_width = 0.2
bar_locations = np.arange(1, 5)
opacity = 1

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

#P = [0.57, 0.62, 0.62, 0.45]
#PT = [0.403, 0.463, 0.44, 0.08]
#PT_lambda = [0.00246, 0.0030, 0.0040, 0.0033]

#P_err = [_, _, 0.032, _]
#PT_err = [0.01644, 0.01582, _, _]
#PT_lambda = [0.0001624, 0.000233, _, _]

#params = np.array([P, PT, PT_lambda])

#ax1.fill_between([1.05, 1.55], [0, 0], [0.7, 0.7], color='grey', alpha=0.3)
#ax1.fill_between([2.05, 2.55], [0, 0], [0.7, 0.7], color='grey', alpha=0.3)
#ax1.fill_between([3.05, 3.55], [0, 0], [0.7, 0.7], color='grey', alpha=0.3)
#ax1.fill_between([4.05, 4.55], [0, 0], [0.7, 0.7], color='grey', alpha=0.3)

for i, res in enumerate(sol):

        params = file_rw.load_object('%s/%s/10wt/forecast_%s_1statepowerlaw_gaussian.pl' % (root_dir, res, res))
        powerlaw_alpha = []
        for j in params.dwell_parameters[0]:
            powerlaw_alpha.append(j)
        
        params = file_rw.load_object('%s/%s/10wt/forecast_%s_1state_powercut_levy.pl' % (root_dir, res, res))
        powercut_alpha = []
        powercut_lambda = []
        for j in params.dwell_parameters[0]:
            powercut_alpha.append(j[0])
            powercut_lambda.append(j[1])

        ax1.bar(bar_locations[i] - bar_width, np.mean(powerlaw_alpha), bar_width, label=names[res], edgecolor='black', color=colors['powerlaw_alpha'], alpha=opacity, yerr=np.std(powerlaw_alpha))
        ax1.bar(bar_locations[i], np.mean(powercut_alpha), bar_width, label=names[res], edgecolor='black', color=colors['powercut_alpha'], alpha=opacity, yerr=np.std(powercut_alpha))
        ax2.bar(bar_locations[i] + bar_width, np.mean(powercut_lambda), bar_width, label=names[res], edgecolor='black', color=colors['powercut_lambda'], alpha=opacity, yerr=np.std(powercut_lambda))
        #ax2.bar(bar_locations[i] + bar_width / 2, H[i], bar_width, hatch='//', edgecolor='black', color=colors[res], alpha=0.7)
	#plt.bar(bar_locations + i*bar_width - bar_width/2, heights, bar_width, label=names[res])

#ax1.text(0.65, 0.65, r'P($\alpha_d$)', fontsize=14)
#ax1.text(1.05, 0.65, r'P($\alpha_d, \lambda$)', fontsize=14)

hatch1 = mpatches.Patch(facecolor=colors['powerlaw_alpha'], label=r'$P(\mathbf{\alpha_d})$', edgecolor='black')
hatch2 = mpatches.Patch(facecolor=colors['powercut_alpha'], label=r'$P_T(\mathbf{\alpha_d}, \lambda)$', edgecolor='black')
hatch3 = mpatches.Patch(facecolor=colors['powercut_lambda'], label=r'$P_T(\alpha_d, \mathbf{\lambda})$', edgecolor='black')

plt.legend(handles=[hatch1, hatch2, hatch3], fontsize=14, loc='upper left')
plt.xticks([1.0, 2.0, 3.0, 4.0], names2)
#for xtick, res in zip(ax1.get_xticklabels(), sol):
#    xtick.set_color(colors[res])

#ax1.set_xlabel('Solute', fontsize=14)
ax1.set_ylabel(r'$\alpha_d$', fontsize=14)
ax2.set_ylabel('$\lambda$', fontsize=14)
ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
ax1.set_ylim(0.0, 1)
#ax2.set_ylim(0.2, 0.45)
plt.tight_layout()
plt.savefig('1mode_AD_dwells.pdf')
plt.show()
