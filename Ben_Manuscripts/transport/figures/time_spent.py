#!/usr/bin/env python

import numpy as np
from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import names
from scipy import stats

def calculate_rdf(res, path, gro='berendsen.gro', traj='PR_nojump.xtc', atoms=None):

	print('Calculating RDF of residue %s' % r)
	if atoms is not None:
		rdf = System('%s/%s' %(path, gro), '%s/%s' %(path, traj), r, 'HII', atoms=atoms)
	else:
		rdf = System('%s/%s' %(path, gro), '%s/%s' %(path, traj), r, 'HII')

	rdf.radial_distribution_function(bins=50, spline=True, npts_spline=10, cut=1.5)

	rdf.bootstrap(200)
	
	file_rw.save_object(rdf, '%s/rdf_%s.pl' % (path, res))

	return rdf

recalculate = False 
residues = ['ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE']

radius = 0.723 # for phenyls
#radius=0.513 # for carboxylates
wt=10

frac = np.zeros([len(residues), 2])

try:
	load = np.load('frac_time_r_%.3f.npz' % radius)
	frac = load['frac']

except FileNotFoundError:
	for i, r in enumerate(residues):

		path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" %(r,wt)

		if recalculate:
			rdf = calculate_rdf(r, path)
		else:
			try:
				rdf = file_rw.load_object('%s/rdf_%s.pl' %(path, r))
			except FileNotFoundError:
				rdf = calculate_rdf(r, path)

		zbox = rdf.t.unitcell_vectors[:, 2, 2].mean()
		rd = rdf.radial_distances

		ntransitions = 0  # number of times solute switches between pore and tail region	
		for res in rd.T:
			inpore = np.zeros_like(res)
			inpore[np.where(res < radius)[0]] = True
			try:
				ntransitions += len(np.argwhere(np.diff(inpore)).squeeze().tolist())  # see forecast_ctrw.py
			except TypeError:
				print(np.argwhere(np.diff(inpore)).squeeze().tolist())
				ntransitions += 1

		#mean = rdf.density.mean(axis=0)
		#V = np.array([zbox * mean[i] * np.pi*(rdf.r[i + 1] ** 2 - rdf.r[i] ** 2) for i in range(len(rdf.r) - 1)])
		#divider = np.argmin(np.abs(rdf.r - radius))
		#f = np.sum(V[:divider]) / np.sum(V)
		f = np.where(rd.flatten() < radius)[0].size / rd.flatten().size
		err = np.sqrt((1 - f) * f / ntransitions)
		print("Fraction of time spent within %.2f nm of pore center by %s : %.3f +\- %.3f" % (radius, r, f, err))
		frac[i] = [f, err]

		np.savez_compressed('frac_time_r_%.3f' % radius, frac=frac)

ordered = np.argsort(frac[:, 0])[::-1]
labels = np.array([names.abbreviation[r] for r in residues])[ordered]
colors = np.array([names.color_dict[r] for r in residues])[ordered]
print(frac)
index = np.arange(len(residues))
fig, ax = plt.subplots()
ax.tick_params(labelsize=14)
ax.set_xticks(index)
ax.set_xticklabels(labels, fontsize=14)
[x.set_color(colors[i]) for i, x in enumerate(plt.gca().get_xticklabels())]
plt.xticks(rotation=90)
ax.bar(index, frac[ordered, 0], yerr=frac[ordered, 1])
ax.set_ylabel('Fraction of time spent in pore region', fontsize=14)
savename = 'frac_time_spent.pdf'
plt.tight_layout()
plt.savefig(savename)
plt.show()


