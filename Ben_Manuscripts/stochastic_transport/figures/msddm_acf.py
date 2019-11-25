#!/usr/bin/env python

from LLC_Membranes.analysis.markov_state_dependent_dynamics import States
from LLC_Membranes.llclib import timeseries, fitting_functions, file_rw
from scipy.stats import levy_stable
import matplotlib.pyplot as plt
import numpy as np

res = 'URE'
directory = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt" % res
ntraj = 24  # number of trajectories to simulate
nboot = 200  # number of bootstrap trials when getting errorbars on MSD
max_k = 15

equil = {'GCL': 2400, 'URE': 2000, 'MET': 7000, 'ACH': 8800}  # frame number, not ns. (multiply ns by 2)
truncate = {'GCL': 0.8, 'URE': 0.8, 'MET': 1.0, 'ACH': 1.0}

traj = '5ms_nojump.xtc'
gro = 'em.gro'

first_frame = equil[res]  # frame at which to start reading trajectory

# probably easier to just re-run these calculations in the appropriate directory. 
# Doesn't matter which dwell/hop is used as they will be re-fit below
states = file_rw.load_object('%s/states.pl' % directory)
nstates = states.count_matrix.shape[0]
ntransitions = 0
for i in range(nstates):
	for j in range(nstates):
		if i != j:
			ntransitions += states.count_matrix[i, j]
print(states.count_matrix)
print(states.count_matrix.sum())
print(ntransitions)
percent_transitions = 100 * (ntransitions / states.count_matrix.sum())
print('Percentage transitions: %.2f' % percent_transitions)
acf = states.transition_autocorrelation()
timeseries.plot_autocorrelation(acf, max_k=max_k, nboot=200, show=False, label='Empirical autocorrelation')
H = timeseries.hurst(acf, nboot=nboot, max_k=max_k).mean()
plt.plot(np.arange(max_k + 1), fitting_functions.hurst_autocovariance(np.arange(max_k + 1), H), '--', color='black', lw=2, label='Fit theoretical autocorrelation')
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('msddm_acf.pdf')
plt.show()
#states.calculate_hurst(plot=True, nboot=nboot, max_k=max_k, state=8)
