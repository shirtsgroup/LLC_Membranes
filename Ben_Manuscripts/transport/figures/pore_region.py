#!/usr/bin/env python

import numpy as np
from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw, stats
import matplotlib.pyplot as plt
import names
import tqdm

"""
Determine boundary of pore region based on average position of some
group of atoms. 
"""

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

pickle = "rdf_HII_CC1C2C3C4C5.pl"
pickle = "rdf_HII_O3O4.pl"

residues = ["ACH", "ACN", "ATO", "BUT", "DMF", "DMP", "DMS", "EAC", "ETH", "GCL", "GLY", "MET", "PCB", "PG", "PR", "RIB", "SOH", "TET", "THF", "URE"] 

nframes = 200
#	d_head_groups = np.zeros([len(residues)*nframes, 50, len(residues)])
d_head_groups = np.zeros([len(residues)*nframes, 50])
wt = 10

rmax = 0.85  # only look for peaks up to this radius. These rdfs generally have two similar-height peaks. We want the r location of the first peak.

maxes = []
for i, r in enumerate(residues):

	path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" %(r,wt)

	hg = file_rw.load_object('%s/%s' % (path, pickle))
	#d_head_groups[:, i] = hg.density.mean(axis=0)
	d_head_groups[i*nframes:(i+1)*nframes, :] = hg.density
	max_ndx = np.argmin(np.abs(hg.r - rmax))
	maxes.append(hg.r[np.argmax(hg.density.mean(axis=0)[:max_ndx])])
	print(maxes)
	
print("Mean Max Head Group Density: %.3f" % np.mean(maxes))

