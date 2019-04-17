#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import names
from LLC_Membranes.analysis import lifetime


# ATO, DMF, DMS, EAC, PCB, THF don't hbond with head groups
residues = ['ACH', 'ACN', 'BUT', 'DMP', 'ETH', 'GCL', 'GLY', 'MET', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'URE']

ci = .95

mean_5wt = np.zeros([len(residues), 2])
ci_5wt = np.zeros([len(residues), 2])
mean_10wt = np.zeros([len(residues), 2])
ci_10wt = np.zeros([len(residues), 2])

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"

try:
	loaded = np.load('hbond_lifetimes_ci_%.2f.npz' % ci)
	mean_5wt = loaded['mean_5wt']
	ci_5wt = loaded['ci_5wt']
	mean_10wt = loaded['mean_10wt']
	ci_10wt = loaded['ci_10wt']

except FileNotFoundError:
	for i, r in enumerate(residues):
		for wt in [10, 5]:
			print(r, '%swt' % wt)
			full_path = '%s/%s/%swt' % (path, r, wt)
			life = lifetime.HydrogenBondLifetime('%s/PR_nojump.xtc' % full_path, '%s/berendsen.gro' % full_path)
			life.set_eligible(r, 'all')
			life.set_eligible('HII', 'all')
			life.identify_hbonds(.35, 30)
			life.identify_unique_hbonds()
			life.create_timeseries()
			life.calculate_lifetimes(ci=ci)
			if wt == 5:
				mean_5wt[i] = life.mean_lifetime
				ci_5wt[i] = life.confidence
			elif wt == 10:
				mean_10wt[i] = life.mean_lifetime
				ci_10wt[i] = life.confidence
			print(life.mean_lifetime, life.confidence)

	np.savez_compressed('hbond_lifetimes_ci_%.2f' % ci, mean_5wt=mean_5wt, ci_5wt=ci_5wt, mean_10wt=mean_10wt, ci_10wt=ci_10wt)

# saved values
#                     ACH   ACN  BUT   DMP   ETH    GCL   GLY    MET   PG    PR    RIB    SOH    TET  URE
#mean_5wt = np.array([5.56, 1.86, 2.10, 1.84, 2.14,  2.07, 3.39,  1.81, 2.69, 2.60, 4.07,  1.62,  3.24, 1.54]) 
#median_5wt = np.array([1,  0.5,  0.5,  0.5,  0.5,   0.5,  0.5,   1,    0.50, 0.50, 0.5,   0.50,  1, 0.5])
#maximum_5wt = np.array([248.5, 47,130, 753,  283.5, 144,  304.50, 68,  457,  225,  538.5, 189.5, 431, 46.50])

#                     ACH   ACN   BUT   DMP   ETH   GCL   GLY   MET   PG    PR    RIB   SOH   TET    URE
#mean_10wt = np.array([4.01, 1.10, 1.20, 1.37, 1.33, 1.34, 2.36, 1.31, 1.95, 1.24, 2.77, 1.21, 2.63,  0.96]) 
#median_10wt = np.array([1,  0.5,  0.50, 0.5,  0.5,  0.5,  0.5,  0.5,  0.50, 0.5,  0.5,  0.50, 0.5,   0.50])
#maximum_10wt = np.array([284,14, 80.50, 105,  44.5, 57.5, 158,  63,   296,  214,  203.5, 82,  379.5, 14.50])

ordered = np.argsort(ci_10wt[:, 0])[::-1]
labels = np.array([names.abbreviation[r] for r in residues])[ordered]
colors = np.array([names.color_dict[r] for r in residues])[ordered]

bar_width = 0.4
opacity = 0.8
index = np.arange(len(residues))
fig, ax = plt.subplots()
ax.tick_params(labelsize=14)
ax.set_xticks(index)
ax.set_xticklabels(labels, fontsize=14)
[x.set_color(colors[i]) for i, x in enumerate(plt.gca().get_xticklabels())]
plt.xticks(rotation=90)
ax.bar(index - bar_width / 2, ci_10wt[ordered, 0], bar_width, alpha=opacity, label='10 wt % water', yerr=ci_10wt[ordered, 1])
ax.bar(index + bar_width / 2, ci_5wt[ordered, 0], bar_width, alpha=opacity, label='5 wt % water', yerr=ci_5wt[ordered, 1])
ax.set_ylabel('95$^{th}$ percentile Hbond Lifetime (ns)', fontsize=14)
savename = 'hbond_lifetime.pdf'
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(savename)
plt.show()
