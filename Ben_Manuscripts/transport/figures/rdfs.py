#!/usr/bin/env python

import numpy as np
from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw, stats
import matplotlib.pyplot as plt
import names
import tqdm

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
simple_alcohols = False
polyols = False 
head_groups = True 
thiol_comparison = False
ketones = False 
nondonors = True
probability = False 

if simple_alcohols:
	residues=["MET", "ETH", "PR", "BUT"]  # simple_alcohol_rdf.pdf 
elif polyols:
	residues=["GCL", "PG", "GLY", "TET", "RIB"]
elif thiol_comparison:
	residues=["SOH", "GCL"]
	#residues=["DMP", "GLY"]
	#residues=["DMS", "ATO"]
elif ketones:
	residues=["ACH", "URE", "ACN", "ATO"]
elif nondonors:
	residues=["THF", "PCB", "EAC", "DMF"]
	#residues=["THF", "DMF"]
else:
	residues=["PG", "GCL"]
	# residues=["DMP", "GLY"]
	#residues = ["GLY", "TET", "RIB"]

#residues = ["BUT", "THF", "PCB", "EAC", "DMF"]
wt=10
maximum = 0
i = 0
v = np.zeros([len(residues), 49])
equil = 200  # chop off first equil frames
for r in residues:

	path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" %(r,wt)

	if recalculate:
		rdf = calculate_rdf(r, path)
	else:
		try:
			rdf = file_rw.load_object('%s/rdf_%s.pl' %(path, r))
		except FileNotFoundError:
			rdf = calculate_rdf(r, path)
	mean = rdf.density[equil:].mean(axis=0)

	if probability:
		rdf.errorbars /= sum(mean)
		mean /= sum(mean)

	new_max = np.amax(mean[np.argwhere(rdf.r > 0.4)])  # really looking for the head group peak
	maximum = max(maximum, new_max)
	plt.plot(rdf.r, mean, label='%s' % names.res_to_name[r])
	plt.fill_between(rdf.r, rdf.errorbars[1, :] + mean, mean - rdf.errorbars[0, :], alpha=0.7)
	#v[i, :] = [mean[i] * np.pi*(rdf.r[i + 1] ** 2 - rdf.r[i] ** 2) for i in range(len(rdf.r) - 1)]
	#print(r, sum(v[i, :np.argmin(np.abs(rdf.r - 0.4)**2)]))
	#plt.plot(rdf.r[:-1], v[i, :])
	i += 1

nboot = 200
if head_groups:

	nframes = 2000
#	d_head_groups = np.zeros([len(residues)*nframes, 50, len(residues)])
	d_head_groups = np.zeros([len(residues)*nframes, 50])

	for i, r in enumerate(residues):
		
		path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" %(r,wt)

		hg = file_rw.load_object('%s/rdf_HII_CC1C2C3C4C5.pl' % path)
		#d_head_groups[:, i] = hg.density.mean(axis=0)
		d_head_groups[i*nframes:(i+1)*nframes, :] = hg.density
	
	boot = np.zeros([nboot, 50])
	for b in tqdm.tqdm(range(nboot)):
		ndx = np.random.choice(np.arange(d_head_groups.shape[0]), size=d_head_groups.shape[0], replace=True)
		boot[b, :] = d_head_groups[ndx, :].mean(axis=0)

	mean = boot.mean(axis=0)
	error = stats.confidence_interval(boot, 95) * (maximum / np.max(mean))
	mean *= (maximum / np.max(mean))
	#std *= (24 / 400)
	#mean *= (24 / 400)
	
	#plt.plot(hg.r, maximum * hg.density.mean(axis=0) / np.max(hg.density.mean(axis=0)), '--')

	# Option 1 
	plt.plot(hg.r, mean, '--', color='black', label='Head Groups')
	#plt.fill_between(hg.r, mean + std, mean - std, alpha=0.6, color='black')
	plt.fill_between(hg.r, mean + error[1, :], mean - error[0, :], alpha=0.6, color='black')

	# Option 2
	#rmax = hg.r[np.argmax(mean)]
	#plt.plot([rmax, rmax], [0, 1.05*mean.max()], '--', color='black')

	# Option 3
	# rmax_ndx = np.argmax(mean)
	# r = hg.r[:rmax_ndx]
	# mean = mean[:rmax_ndx]
	# std = std[:rmax_ndx]
	# plt.plot(r, mean, '--', color='black')
	#plt.fill_between(r, mean + std, mean - std, alpha=0.7, color='black')



#r = residues[0]
#path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" % (r, wt)

#try:
#	rdf = file_rw.load_object('%s/rdf_HII.pl' %path)

#except FileNotFoundError:

#	rdf = calculate_rdf(r, path, atoms=['C', 'C1', 'C2', 'C3', 'C4', 'C5'])

#normalization = 24 / 400
#plt.plot(rdf.r, maximum * rdf.density.mean(axis=0) / np.amax(rdf.density.mean(axis=0)), '--', color='black')
#plt.plot(rdf.r, normalization * rdf.density.mean(axis=0), '--', color='black')

plt.ylabel('Density (count / nm$^3$)', fontsize=14)
plt.xlabel('Distance from pore center (nm)', fontsize=14)
#plt.ylim(-0.05, 1.3) # for diols only
#plt.ylim(-0.015, 0.45) # for DMSO and acetone thiol comparison
# plt.ylim(-0.015, 0.5) # for mercaptoethanol and ethylene glycol comparison
plt.ylim(-0.015, 0.55) # for nondonors
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=13, loc=1, ncol=2, columnspacing=0.5)
plt.tight_layout()
if simple_alcohols:
	plt.savefig('simple_alcohol_rdf.pdf')
elif polyols:
	plt.savefig('polyols_rdf.pdf')
elif thiol_comparison:
	plt.savefig('thiol_comparison_%s.pdf' % residues[0])
elif ketones:
	plt.savefig('ketone_rdf.pdf')
elif nondonors:
	plt.savefig('nondonors_rdf.pdf')
plt.show()

