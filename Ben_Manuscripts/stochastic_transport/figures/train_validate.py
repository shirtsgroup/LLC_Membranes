#!/usr/bin/env python

from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
from LLC_Membranes.timeseries.msd import Diffusivity
from LLC_Membranes.llclib import file_rw
from LLC_Membranes.timeseries.ctrwsim import CTRW
import matplotlib.pyplot as plt
import numpy as np
import argparse

def initialize():

	parser = argparse.ArgumentParser()

	parser.add_argument('-r', '--residue')
	parser.add_argument('-m', '--nmodes', default=1, type=int)
	parser.add_argument('-p', '--portion', default='front', help='front or back')
	parser.add_argument('-msd', '--msd', action="store_true")

	return parser

def md_diffusivity(front=True):

	MD_MSD = Diffusivity('%s/%s' % (directory,traj), '%s/%s' % (directory, gro), 'z', begin=first_frame, end=last_frame, residue=res)

	MD_MSD.calculate()
	MD_MSD.bootstrap(200, fit_line=False)

	# reduce object size and save
	MD_MSD.t = None
	MD_MSD.com = None
	file_rw.save_object(MD_MSD, '%s_msd_train_%s.pl' % (res, suffix))

	return MD_MSD

args = initialize().parse_args()

if args.residue is None:
	import sys
	sys.exit('Please tell me which residue to look at. For example:\n./train_test.py -r URE')

res = args.residue
directory = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt" % res
fracshow = 0.4  # fraction of MD MSD to plot
recalculate_msd = args.msd 
recalculate_walks = False 

front = False
if args.portion.lower() == 'front':
	front = True  # if True, use the parameters trained on the first half of the dataset. Otherwise use second half.
endframe = 1000

nmodes = args.nmodes 
ntraj = 1000  # number of trajectories to simulate
nboot = 200  # number of bootstrap trials when getting errorbars on MSD
padding = 10  # higher number gives better resolution to CTRW trajectories
nt = 8

equil = {'GCL': 2400, 'URE': 2000, 'MET': 7000, 'ACH': 4000}  # frame number, not ns. (multiply ns by 2)
nframes = {'URE': 10131, 'GCL': 10165, 'MET': 11808, 'ACH': 12468}
traj = '5ms_nojump.xtc'
gro = 'em.gro'

tot_frames = nframes[res] - equil[res]
middle_frame = equil[res] + tot_frames // 2

if front: 
	suffix = 'front'  # parameters to use
	# msd-specific stuff
	msd_suffix = 'back'
	first_frame = middle_frame
	last_frame = -1
	last_frame = nframes[res]

else:
	suffix = 'back'
	msd_suffix = 'front'
	first_frame = equil[res]
	last_frame = middle_frame

print(first_frame, last_frame)

if recalculate_msd:

	MD_MSD = md_diffusivity()
else:
	try:
		MD_MSD = file_rw.load_object('%s_msd_train_%s.pl' % (res, msd_suffix))

	except FileNotFoundError:

		MD_MSD = md_diffusivity()

#endshow = int(MD_MSD.nT * fracshow)
endshow = endframe
dt = MD_MSD.dt
plt.plot(np.arange(endshow)*dt, MD_MSD.MSD_average[:endshow], color='black', lw=2)
plt.fill_between(np.arange(endshow)*dt, MD_MSD.MSD_average[:endshow] + MD_MSD.limits[0, :endshow], MD_MSD.MSD_average[:endshow] - MD_MSD.limits[1, :endshow], alpha=0.3, color='black')

last = MD_MSD.MSD_average[endshow]
print('MD MSD: %.2f [%.2f, %.2f]' % (last, last - MD_MSD.limits[1, endshow], last + MD_MSD.limits[0, endshow]))
labels = ['MD']

nsteps = MD_MSD.nT  # match the number of frames i
#endframe = int(fracshow * padding * nsteps)
endframe *= padding

# probably easier to just re-run these calculations in the appropriate directory. 
# Doesn't matter which dwell/hop is used as they will be re-fit below
sys = file_rw.load_object('%s/forecast_%s_%dstate_train_%s.pl' % (directory, res, nmodes, suffix))
dwell_dist = ['Power Law', 'Power Law Exponential Cutoff', 'Power Law', 'Power Law Exponential Cutoff']
hop_dist = ['Gaussian', 'Gaussian', 'Levy', 'Levy']

#dwell_dist = ['Power Law Exponential Cutoff', 'Power Law Exponential Cutoff']
#hop_dist = ['Gaussian', 'Levy']

#dwell_dist = ['Power Law Exponential Cutoff']
#dwell_dist = ['Power Law']
#hop_dist = ['Levy']

# names here will appear in legend
abbreviations = {'Power Law Gaussian': 'sFBM', 'Power Law Exponential Cutoff Gaussian': 'sFBMcut', 'Power Law Levy': 'sFLM', 'Power Law Exponential Cutoff Levy': 'sFLMcut'}

for i, dists in enumerate(zip(dwell_dist, hop_dist)):

	dwell = dists[0]
	hop = dists[1]

	abbrev = abbreviations[dwell + ' ' + hop]
	print(abbrev)

	savename = '%s_%s_%dmode_train_%s.pl' % (res, abbrev, nmodes, suffix)
	#savename = '%s_%s.pl' % (res, abbrev)

	if recalculate_walks:
		sys.fit_distributions(nbins=50, nboot=nboot, plot=False, show=False, save=False, dwell_distribution=dwell, hop_distribution=hop)
	else:
		sys, random_walks = file_rw.load_object(savename)
	
	#print('%s Hop Distribution Parameters:' % abbrev)
	if hop == 'Gaussian':
		motion = 'fbm'
		for m in range(nmodes):
			mean_sigma = np.mean([p[1] for p in sys.hop_parameters[m]])
			mean_mu = np.mean([p[0] for p in sys.hop_parameters[m]])
			#print('Mode %d' % (m + 1))
			#print('mu = %.2f\nsigma = %.2f' % (mean_mu, mean_sigma))
	elif hop == 'Levy':
		motion = 'flm'
		for m in range(nmodes):
			mean_alpha = np.mean([p[0] for p in sys.hop_parameters[m]])
			mean_sigma = np.mean([p[2] for p in sys.hop_parameters[m]])
			mean_mu = np.mean([p[1] for p in sys.hop_parameters[m]])
			#print('Mode %d' % (m + 1))
			#print('alpha = %.3f\nsigma = %.2f\nmu = %.2f' % (mean_alpha, mean_sigma, mean_mu))
	else:
		import sys
		sys.exit('Type of motion undefined')

	#print('%s Dwell Time Distribution Parameters:' % abbrev)
	if dwell == 'Power Law':
		for m in range(nmodes):
			mean_alpha = np.mean(sys.dwell_parameters[m])
			#print('Mode %d' % (m + 1))
			#print('alpha=%.4f' % mean_alpha)
	elif dwell == 'Power Law Exponential Cutoff':
		for m in range(nmodes):
			mean_alpha = np.mean([p[0] for p in sys.dwell_parameters[m]])
			mean_lambda = np.mean([p[1] for p in sys.dwell_parameters[m]])
			#print('Mode %d' % (m + 1))
			#print('alpha=%.4f\nlambda=%.4f' % (mean_alpha, mean_lambda))

	print('\nHurst parameter: %.2f' % np.mean(sys.hurst_distribution))
	if recalculate_walks:
		if nmodes > 1:
			sys.determine_transition_matrix()
		random_walks = CTRW(nsteps, ntraj, nmodes=nmodes, dt=dt, hop_dist=motion, dwell_dist=dwell, transition_count_matrix=sys.count_matrix if sys.nmodes > 1 else None)
		random_walks.generate_trajectories(fixed_time=True, distributions=(sys.dwell_parameters, sys.hop_parameters, sys.hurst_distribution), discrete=True, ll=sys.dwell_lower_limit, max_hop=2*sys.max_hop, nt=nt)

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
#plt.savefig('toc_msd.png')
plt.savefig('%dmode_msd_comparison_%s_train_%s.pdf' % (nmodes,res,suffix))
plt.show()
