#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.timeseries.forecast_ctrw import System
from LLC_Membranes.llclib import file_rw
import names

residues = ["GCL", "SOH"]
wt = 10
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
colors = ['blue', 'red']
opacity = 1 
nbins = 25
lw = 2

fig, ax = plt.subplots(1, 2, figsize=(10, 5))

for j, r in enumerate(residues):
	obj = file_rw.load_object('%s/%s/%swt/forecast_%s.pl' % (path, r, wt, r))
	hops = []
	for i in obj.hop_lengths:
		hops += i
	print(max(hops))
	if j == 0:
		hop_hist, edges = np.histogram(hops, density=True, bins=nbins)
		bounds = [edges[0], edges[-1]]
	else:
		hop_hist, edges = np.histogram(hops, density=True, bins=np.linspace(bounds[0], bounds[1], nbins + 1)) 

	hop_outline = np.zeros([len(hop_hist)*2 + 2, 2])
	hop_outline[::2, 0] = edges
	hop_outline[1::2, 0] = edges
	hop_outline[1:-1:2, 1] = hop_hist
	hop_outline[2:-1:2, 1] = hop_hist

	if j == 0:
		dwell_hist, edges = np.histogram(obj.dwell_times, density=True, bins=nbins)
		bounds_power = [edges[0], edges[-1]]
	else:
		dwell_hist, edges = np.histogram(obj.dwell_times, density=True, bins=np.linspace(bounds_power[0], bounds_power[1], nbins + 1))

	dwell_outline = np.zeros([len(dwell_hist)*2 + 2, 2])
	dwell_outline[::2, 0] = edges
	dwell_outline[1::2, 0] = edges
	dwell_outline[1:-1:2, 1] = dwell_hist
	dwell_outline[2:-1:2, 1] = dwell_hist

	ax[0].plot(hop_outline[:, 0], hop_outline[:, 1], color=colors[j], alpha=opacity, linewidth=lw)
	ax[1].plot(dwell_outline[:, 0], dwell_outline[:, 1], color=colors[j], alpha=opacity, label=names.res_to_name[r], linewidth=lw)

ax[0].tick_params(labelsize=14)
ax[1].tick_params(labelsize=14)
ax[1].legend(fontsize=14)
ax[0].set_ylabel('Frequency', fontsize=14)
ax[0].set_xlabel('Hop Length (nm)', fontsize=14)
ax[1].set_xlabel('Dwell Time (ns)', fontsize=14)
plt.tight_layout()
plt.savefig('dwell_hop_%s.pdf' % '_'.join(residues))
plt.show()
