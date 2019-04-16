#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import names
from LLC_Membranes.analysis.hop_location import HopLocation
from LLC_Membranes.llclib import file_rw
import argparse

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
wt = 10
residues = ['ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE']

parser = argparse.ArgumentParser()

parser.add_argument('-c', '--cut', default=0.5, type=float, help='Cut-off distance between pores and tails')

args = parser.parse_args()

cut = args.cut

print('Cut-off = %s' % cut)

try:

	loaded = np.load('hop_lengths_cut_%s.npz' % cut)
	inpore = loaded['inpore']
	outpore = loaded['outpore']

except FileNotFoundError:

	inpore = np.zeros([len(residues), 2]) # [nresidues, (mean, std)]
	outpore = np.zeros([len(residues), 2])

	for i, r in enumerate(residues):
		print(r)

		hops = file_rw.load_object('%s/%s/%swt/hop_locations.pl' % (path, r, wt))
		x = hops.regional_hop_lengths(cut)
		inpore[i, 0] = x[0]
		inpore[i, 1] = x[1]
		outpore[i, 0] = x[2]
		outpore[i, 1] = x[3]

	np.savez_compressed('hop_lengths_cut_%s' % cut, inpore=inpore, outpore=outpore)

# saved values
#inpore = np.array([0.25, 0.23, 0.21, 0.18, 0.16, 0.16, 0.24, 0.19, 0.28, 0.30, 0.20, 0.39, 0.18, 0.22, 0.23, 0.15, 0.27, 0.21, 0.22, 0.26])
#outpore = np.array([0.17, 0.18, 0.15, 0.11, 0.14, 0.13, 0.15, 0.13, 0.19, 0.21, 0.15, 0.30, 0.10, 0.16, 0.17, 0.09, 0.18, 0.14, 0.15, 0.17])

percent = (inpore[:, 0] - outpore[:, 0]) / outpore[:, 0]
print('On average, solute hop lengths are %.2f %% larger in the pore region' % (100*percent.mean()))

ordered = np.argsort(inpore[:, 0])[::-1]
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
ax.bar(index - bar_width / 2, inpore[ordered, 0], bar_width, alpha=opacity, label='Hops in pore region', yerr=inpore[ordered, 1])
ax.bar(index + bar_width / 2, outpore[ordered, 0], bar_width, alpha=opacity, label='Hops outside pore region', yerr=outpore[ordered, 1])
ax.set_ylabel('Mean hop length (nm)', fontsize=14)
savename = 'hop_length.pdf'
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig(savename)
plt.show()
