#!/usr/bin/env python

from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
from LLC_Membranes.analysis.msd import Diffusivity
from LLC_Membranes.llclib import file_rw
from LLC_Membranes.timeseries.ctrwsim import CTRW
import matplotlib.pyplot as plt
import numpy as np

def md_diffusivity():

        MD_MSD = Diffusivity('%s/%s' % (directory,traj), '%s/%s' % (directory, gro), 'z', begin=first_frame, residue=res)
        MD_MSD.calculate()
        MD_MSD.bootstrap(nboot, fit_line=False)
        MD_MSD.t = None
        MD_MSD.com = None

	# reduce object size and save
        MD_MSD.t = None
        MD_MSD.com = None
        file_rw.save_object(MD_MSD, '%s_msd.pl' % res)

        return MD_MSD


res = 'URE'
directory = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt" % res
fracshow = 0.4  # fraction of MD MSD to plot
recalculate_msd = False 
recalculate_walks = True
nmodes = 1 
ntraj = 1000  # number of trajectories to simulate
nboot = 200  # number of bootstrap trials when getting errorbars on MSD
padding = 10  # higher number gives better resolution to CTRW trajectories

equil = {'GCL': 2400, 'URE': 2000, 'MET': 7000, 'ACH': 8800}  # frame number, not ns. (multiply ns by 2)

traj = '5ms_nojump.xtc'
gro = 'em.gro'

first_frame = equil[res]  # frame at which to start reading trajectory

if recalculate_msd:

	MD_MSD = md_diffusivity()
else:
	try:
		MD_MSD = file_rw.load_object('%s_msd.pl' % res)

	except FileNotFoundError:

		MD_MSD = md_diffusivity()
	
endshow = int(MD_MSD.nT * fracshow)
dt = MD_MSD.dt
plt.plot(np.arange(endshow)*dt, MD_MSD.MSD_average[:endshow], color='black', lw=2)
plt.fill_between(np.arange(endshow)*dt, MD_MSD.MSD_average[:endshow] + MD_MSD.limits[0, :endshow], MD_MSD.MSD_average[:endshow] - MD_MSD.limits[1, :endshow], alpha=0.3, color='black')
labels = ['MD']

nsteps = MD_MSD.nT  # match the number of frames 
endframe = int(fracshow * padding * nsteps)

# probably easier to just re-run these calculations in the appropriate directory. 
# Doesn't matter which dwell/hop is used as they will be re-fit below
sys = file_rw.load_object('%s/forecast_%s_%dstate.pl' % (directory, res, nmodes))
dwell_dist = ['Power Law', 'Power Law Exponential Cutoff', 'Power Law', 'Power Law Exponential Cutoff']
hop_dist = ['Gaussian', 'Gaussian', 'Levy', 'Levy']

# names here will appear in legend
abbreviations = {'Power Law Gaussian': 'sFBM', 'Power Law Exponential Cutoff Gaussian': 'sFBMcut', 'Power Law Levy': 'sFLM', 'Power Law Exponential Cutoff Levy': 'sFLMcut'}

for i, dists in enumerate(zip(dwell_dist, hop_dist)):

	dwell = dists[0]
	hop = dists[1]

	abbrev = abbreviations[dwell + ' ' + hop]
	savename = '%s_%s_%dmode.pl' % (res, abbrev, nmodes)
	#savename = '%s_%s.pl' % (res, abbrev)

	if recalculate_walks:
		sys.fit_distributions(nbins=50, nboot=nboot, plot=False, show=False, save=False, dwell_distribution=dwell, hop_distribution=hop)
	else:
		sys, random_walks = file_rw.load_object(savename)
	
	print('%s Hop Distribution Parameters:' % abbrev)
	if hop == 'Gaussian':
		motion = 'fbm'
		for m in range(nmodes):
			mean_sigma = np.mean([p[1] for p in sys.hop_parameters[m]])
			mean_mu = np.mean([p[0] for p in sys.hop_parameters[m]])
			print('Mode %d' % (m + 1))
			print('mu = %.2f\nsigma = %.2f' % (mean_mu, mean_sigma))
	elif hop == 'Levy':
		motion = 'flm'
		for m in range(nmodes):
			mean_alpha = np.mean([p[0] for p in sys.hop_parameters[m]])
			mean_sigma = np.mean([p[2] for p in sys.hop_parameters[m]])
			mean_mu = np.mean([p[1] for p in sys.hop_parameters[m]])
			print('Mode %d' % (m + 1))
			print('alpha = %.3f\nsigma = %.2f\nmu = %.2f' % (mean_alpha, mean_sigma, mean_mu))
	else:
		import sys
		sys.exit('Type of motion undefined')

	print('%s Dwell Time Distribution Parameters:' % abbrev)
	if dwell == 'Power Law':
		for m in range(nmodes):
			mean_alpha = np.mean(sys.dwell_parameters[m])
			print('Mode %d' % (m + 1))
			print('alpha=%.4f' % mean_alpha)
	elif dwell == 'Power Law Exponential Cutoff':
		for m in range(nmodes):
			mean_alpha = np.mean([p[0] for p in sys.dwell_parameters[m]])
			mean_lambda = np.mean([p[1] for p in sys.dwell_parameters[m]])
			print('Mode %d' % (m + 1))
			print('alpha=%.4f\nlambda=%.4f' % (mean_alpha, mean_lambda))

	print('Hurst parameter: %.2f' % np.mean(sys.hurst_distribution))

	if recalculate_walks:
		if nmodes > 1:
			sys.determine_transition_matrix()
			print(sys.count_matrix)
		random_walks = CTRW(nsteps, ntraj, nmodes=nmodes, dt=dt, hop_dist=motion, dwell_dist=dwell, transition_count_matrix=sys.count_matrix if sys.nmodes > 1 else None)
		random_walks.generate_trajectories(fixed_time=True, distributions=(sys.dwell_parameters, sys.hop_parameters, sys.hurst_distribution), discrete=True, ll=sys.dwell_lower_limit, max_hop=2*sys.max_hop)

		random_walks.calculate_msd(ensemble=False)
		random_walks.bootstrap_msd(fit_linear=False)
		random_walks.trajectories = None
		random_walks.trajectory_hops = None
		random_walks.z_interpolated = None
		file_rw.save_object((sys, random_walks), savename)

	random_walks.plot_msd(show=False, end_frame=endframe, newfig=False)
		
	labels.append(abbrev)

plt.legend(labels, loc=0, fontsize=14)
plt.tight_layout()
plt.savefig('%dmode_msd_comparison_%s.pdf' % (nmodes,res))
plt.show()
