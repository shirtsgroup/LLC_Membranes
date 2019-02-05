#!/usr/bin/env python

import numpy as np
from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import names

def calculate_rdf(res, path, gro='berendsen.gro', traj='PR_nojump.xtc', atoms=None):

	print('Calculating RDF of residue %s' % r)
	rdf = System('%s/%s' %(path, gro), '%s/%s' %(path, traj), r, 'HII')

	if atoms is not None:
		rdf.radial_distribution_function(bins=50, spline=True, npts_spline=20, cut=1.5, atoms=atoms)
	else:
		rdf.radial_distribution_function(bins=50, spline=True, npts_spline=20, cut=1.5)

	rdf.bootstrap(200)
	
	file_rw.save_object(rdf, '%s/rdf_%s.pl' % (path, res))

	return rdf

residues=["BUT", "ETH", "PR", "MET"]  # simple_alcohol_rdf.pdf 
wt=10
maximum = 0
for r in residues:

	path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" %(r,wt)
	try:
		rdf = file_rw.load_object('%s/rdf_%s.pl' %(path, r))
	
	except FileNotFoundError:

		rdf = calculate_rdf(r, path)
		

	mean = rdf.density.mean(axis=0)
	maximum = max(maximum, np.amax(mean))
	plt.plot(rdf.r, mean, label='%s' % names.res_to_name[r])
	plt.fill_between(rdf.r, rdf.errorbars[1, :] + mean, mean - rdf.errorbars[0, :], alpha=0.7)

r = residues[0]
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" % (r, wt)

try:
	rdf = file_rw.load_object('%s/rdf_HII.pl' %path)

except FileNotFoundError:

	rdf = calculate_rdf(r, path, atoms=['C', 'C1', 'C2', 'C3', 'C4', 'C5'])

normalization = 24 / 400
plt.plot(rdf.r, maximum * rdf.density.mean(axis=0) / np.amax(rdf.density.mean(axis=0)), '--', color='black')
#plt.plot(rdf.r, normalization * rdf.density.mean(axis=0), '--', color='black')

plt.ylabel('Density (count / nm$^3$)', fontsize=14)
plt.xlabel('Distance from pore center (nm)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('simple_alcohol_rdf.pdf')
plt.show()

