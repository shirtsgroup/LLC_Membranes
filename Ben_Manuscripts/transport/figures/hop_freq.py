#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.analysis import hop_location
import names
from scipy import stats

residues = ['ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE']

ci = .95
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
wt = 10

hopfreq = np.zeros([len(residues), 2])

try:
	loaded = np.load('hopfreq_ci_%.2f.npz' % ci)
	hopfreq = loaded['hopfreq']

except FileNotFoundError:
	for i, r in enumerate(residues):
		print(r)
		full_path = '%s/%s/%swt' % (path, r, wt)
		hops = hop_location.HopLocation('%s/PR_nojump.xtc' % full_path, '%s/berendsen.gro' % full_path, r)
		hops.calculate_solute_partition(r=1.5, spline=True, membrane_residue='HII', write_tcl=False)
		hops.hops_and_dwells(penalty=0.25, nframes_dwell=10, locations=True)
		hopfreq[i] = hops.hop_frequency()

	np.savez_compressed('hopfreq_ci_%.2f' % ci, hopfreq=hopfreq)

ordered = np.argsort(hopfreq[:, 0])[::-1]
labels = np.array([names.abbreviation[r] for r in residues])[ordered]
colors = np.array([names.color_dict[r] for r in residues])[ordered]

# Use poisson pmf to generate errorbars instead
nhops = hopfreq[ordered, 0] * 1000 # total hops over 1000 ns
err = np.zeros([2, nhops.size])
confidence = 68.27
lower_confidence = (100 - confidence) / 2 
upper_confidence = 100 - lower_confidence
lower_confidence /= 100
upper_confidence /= 100

for i, n in enumerate(nhops):
	k = np.arange(0, 1000)
	cdf = stats.poisson.cdf(k, n)
	bottom = np.argmin(np.abs(cdf - lower_confidence))
	top = np.argmin(np.abs(cdf - upper_confidence))
	err[:, i] = [k[bottom], k[top]]

err /= 1000  # convert to hops per 1 ns
err[0, :] = hopfreq[ordered, 0] - err[0, :]
err[1, :] -= hopfreq[ordered, 0]

bar_width = 0.8
opacity = 0.8
index = np.arange(len(residues))
fig, ax = plt.subplots()
ax.tick_params(labelsize=14)
ax.set_xticks(index)
ax.set_xticklabels(labels, fontsize=14)
[x.set_color(colors[i]) for i, x in enumerate(plt.gca().get_xticklabels())]
plt.xticks(rotation=90)
#ax.bar(index, hopfreq[ordered, 0], bar_width, alpha=opacity, yerr=hopfreq[ordered, 1])
ax.bar(index, hopfreq[ordered, 0], bar_width, alpha=opacity, yerr=err)
ax.set_ylabel('Hop Frequency (ns$^{-1}$)', fontsize=14)
savename = 'hopfreq_total.pdf'
plt.tight_layout()
plt.savefig(savename)
plt.show()

