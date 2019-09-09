#!/usr/bin/env python

from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
from LLC_Membranes.analysis.msd import Diffusivity
from LLC_Membranes.llclib import file_rw
from LLC_Membranes.timeseries.ctrwsim import CTRW
import matplotlib.pyplot as plt
import numpy as np

res = 'URE'
directory = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt" % res
MDMSD = True 

equil = {'GCL': 1200, 'URE': 1000}

traj = '5ms_nojump.xtc'
gro = 'em.gro'

first_frame = equil[res]

if MDMSD:

	MD_MSD = Diffusivity('%s/%s' % (directory,traj), '%s/%s' % (directory, gro), 'z', begin=first_frame, residue=res)
	MD_MSD.calculate()
	MD_MSD.bootstrap(200, fit_line=False)
	endshow = 4000
	plt.plot(np.arange(endshow)*0.5, MD_MSD.MSD_average[:endshow], color='black', lw=2)
	plt.fill_between(np.arange(endshow)*0.5, MD_MSD.MSD_average[:endshow] + MD_MSD.limits[0, :endshow], MD_MSD.MSD_average[:endshow] - MD_MSD.limits[1, :endshow], alpha=0.3, color='black')
	labels = ['MD']
else:
	labels = []

nmodes = 2
nsteps = 10000  # 5000 ns with a 0.5 ns timestep
ntraj = 1000
dt = 0.5  # ns
padding = 10
endframe = int(0.4 * padding * nsteps)
nboot = 200

sys = file_rw.load_object('%s/forecast_%s_%dstate.pl' % (directory, res, nmodes))

# names here will appear in legend
dwell_dist = ['Power Law', 'Power Law Exponential Cutoff', 'Power Law', 'Power Law Exponential Cutoff']
hop_dist = ['Gaussian', 'Gaussian', 'Levy', 'Levy']
#dwell_dist = ['Power Law Exponential Cutoff']
#hop_dist = ['Levy']

for i, dists in enumerate(zip(dwell_dist, hop_dist)):

	dwell = dists[0]
	hop = dists[1]

	sys.fit_distributions(nbins=50, nboot=nboot, plot=False, show=False, save=False, dwell_distribution=dwell, hop_distribution=hop)
	
	if hop == 'Gaussian':
		motion = 'fbm'
	elif hop == 'Levy':
		motion = 'flm'
	else:
		import sys
		sys.exit('Type of motion undefined')

	print(len(sys.hurst_distribution))

	random_walks = CTRW(nsteps, ntraj, nmodes=nmodes, dt=dt, hop_dist=motion, dwell_dist=dwell, transition_matrix=sys.transition_matrix if sys.nmodes > 1 else None)
	random_walks.generate_trajectories(fixed_time=True, distributions=(sys.dwell_parameters, sys.hop_parameters, sys.hurst_distribution), discrete=True, ll=sys.dwell_lower_limit, max_hop=sys.max_hop)

	random_walks.calculate_msd(ensemble=False)
	random_walks.bootstrap_msd(fit_linear=False)
	random_walks.plot_msd(show=False, end_frame=endframe, newfig=False)

	labels.append(dwell + ' + ' + hop)
#labels.append('sFBM')

plt.legend(labels, loc=0, fontsize=14)
plt.tight_layout()
plt.savefig('%dmode_msd_comparison_%s.pdf' % (nmodes,res))
#plt.savefig('sflm_msd_comparison_URE.pdf')
plt.show()
