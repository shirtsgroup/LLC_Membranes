#!/usr/bin/env python

from LLC_Membranes.analysis.markov_state_dependent_dynamics import States, Chain
from LLC_Membranes.timeseries.msd import Diffusivity
from LLC_Membranes.llclib import file_rw
from scipy.stats import levy_stable
import matplotlib.pyplot as plt
import numpy as np
import argparse

def initialize():

	parser = argparse.ArgumentParser(description='Model a continuous time random walk')

	parser.add_argument('-r', '--residue', default='URE')
	parser.add_argument('-n', '--ntraj', default=1000, type=int)

	return parser


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

args = initialize().parse_args()

res = args.residue
directory = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt" % res
fracshow = 0.4  # fraction of MD MSD to plot
endframe = 2000
recalculate_msd = False 
recalculate_walks = False 
extreme_trapping = False 
ntraj = args.ntraj # number of trajectories to simulate
nboot = 200  # number of bootstrap trials when getting errorbars on MSD
nt = 8 # number of threads
load = True 

equil = {'GCL': 2400, 'URE': 2000, 'MET': 7000, 'ACH': 4000}  # frame number, not ns. (multiply ns by 2)
truncate = {'GCL': 1.0, 'URE': 1.0, 'MET': 1.0, 'ACH': 1.0}
m = 256
Mlowerbound = 6000
#M = {'GCL': 8619, 'URE': 8253, 'MET': 11577, 'ACH': 7917}

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
	
#endshow = int(MD_MSD.nT * fracshow)
fracshow = endframe / MD_MSD.nT 
endshow = endframe
dt = MD_MSD.dt
plt.plot(np.arange(endshow)*dt, MD_MSD.MSD_average[:endshow], color='black', lw=2)
plt.fill_between(np.arange(endshow)*dt, MD_MSD.MSD_average[:endshow] + MD_MSD.limits[0, :endshow], MD_MSD.MSD_average[:endshow] - MD_MSD.limits[1, :endshow], alpha=0.3, color='black')
labels = ['MD']

print(MD_MSD.MSD_average[endshow])

nsteps = MD_MSD.nT  # match the number of frames 

if not load:
	# probably easier to just re-run these calculations in the appropriate directory. 
	# Doesn't matter which dwell/hop is used as they will be re-fit below
	states = file_rw.load_object('%s/states.pl' % directory)
	#print(states.hurst)
	#exit()

	if extreme_trapping:
		states.hurst[:, :] = 0
	print(states.hurst.mean(axis=1))
	print(states.fit_params)
	print(ntraj)

	chains = Chain(states.count_matrix, states.fit_params, hurst_parameters=states.hurst, emission_function=levy_stable)
	#chains = Chain(states.count_matrix, states.fit_params, hurst_parameters=None, emission_function=levy_stable)
	chains.generate_realizations(ntraj, nsteps, bound=truncate[res], m=m, Mlowerbound=Mlowerbound, nt=nt)
else:
	chains = file_rw.load_object('%s_msddm_chains1.pl' % res)

chains.calculate_msd()
chains.plot_msd(cutoff=fracshow, label='MSDDM', overlay=True, show=False)
labels.append('MSDDM')
print(chains.msds.mean(axis=1)[int(fracshow*nsteps)])
plt.legend(labels, loc=0, fontsize=14)
plt.tight_layout()
savename = '%s_msddm' % res
if extreme_trapping:
	savename += '_zeroH'
plt.savefig('%s.pdf' % savename)
plt.show()
