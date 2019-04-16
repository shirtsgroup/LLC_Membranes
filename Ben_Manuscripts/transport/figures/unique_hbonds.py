#!/usr/bin/env python

import numpy as np
from LLC_Membranes.analysis.hbonds import System
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import names

nsolute = 24
residues = ['ACH', 'ACN', 'BUT', 'DMP', 'ETH', 'GCL', 'GLY', 'MET', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'URE']
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"

regenerate = False


try:
	loaded = np.load('unique_hbonds.npz')
	n_5wt = loaded['n_5wt']
	n_10wt = loaded['n_10wt']

except FileNotFoundError:

	n_5wt = np.zeros([len(residues), 2])
	n_10wt = np.zeros([len(residues), 2])

	for wt in [5, 10]:

		for i, r in enumerate(residues):
	
			loc = '%s/%s/%swt' % (path, r, wt)
			print(loc)
			sys = file_rw.load_object('%s/hbonds.pl' %loc)
	
			n_unique_donors = []

			res_no, nres = sys.number_residues(r)
	
			for frame in range(sys.t.n_frames):
				res_numbers = [res_no[i] for i in sys.hbonds[frame][0]]
				n_unique_donors.append(len(np.unique(res_numbers)))

			# bootstrap
			n_unique_donors = np.array(n_unique_donors)
			nboot = 200
			n = []
			for b in range(nboot):
				choice = np.random.choice(n_unique_donors, size=n_unique_donors.size, replace=True)
				n.append(100*np.mean(n_unique_donors[choice]) / nsolute)
			print(n)
			exit()

			if wt == 5:
				n_5wt[i] = [100*np.mean(n_unique_donors) / nsolute, np.std(n)]
			elif wt == 10:
				n_10wt[i] = [100*np.mean(n_unique_donors) / nsolute, np.std(n)]

	np.savez_compressed('unique_hbonds.npz', n_5wt=n_5wt, n_10wt=n_10wt)

#else:
	#                  ACH   ACN   BUT   DMP   ETH   GCL   GLY   MET   PG    PR    RIB   SOH TET   URE
#	n_10wt = np.array([52.7, 4.02, 38.9, 48.5, 42.2, 64.4, 78.8, 43.2, 67.5, 39.5, 84.6, 48, 79.6, 4.66])
#	n_5wt = np.array([37.5, 4.26, 33.9, 51.9, 37.3, 52.9, 67.4, 36.2, 61.0, 43.1, 80.3, 40.6, 67.9, 8.86])

ordered = np.argsort(n_10wt[:, 0])[::-1]
labels = np.array([names.abbreviation[r] for r in residues])[ordered]
colors = np.array([names.color_dict[r] for r in residues])[ordered]
fig, ax = plt.subplots()
index = np.arange(n_10wt.shape[0])
bar_width = 0.4
opacity = 0.8
ax.bar(index - bar_width / 2, n_10wt[ordered, 0], bar_width, alpha=opacity, label='10 wt % water', yerr=n_10wt[ordered, 1])
ax.bar(index + bar_width / 2, n_5wt[ordered, 0], bar_width, alpha=opacity, label='5 wt % water', yerr=n_5wt[ordered, 1])
ax.tick_params(labelsize=14)
ax.set_xticks(index)
ax.set_xticklabels(labels, fontsize=14)
ax.set_ylabel('% solutes donating hbonds to monomer', fontsize=14)
[x.set_color(colors[i]) for i, x in enumerate(plt.gca().get_xticklabels())]
plt.xticks(rotation=90)
plt.legend(fontsize=14)
plt.tight_layout()
savename = 'nhbonds_all.pdf'
#plt.savefig(savename)
plt.show()

