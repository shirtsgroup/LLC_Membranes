#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.timeseries.ctrwsim import CTRW
from LLC_Membranes.llclib import stats

# some typical parameters characteristic of solutes in my simulations
sigma = 0.3 # Hop length standard deviation (nm)
alpha = 0.5 # anomalous exponent
hurst = 0.35 # hurst parameter

# Define trajectories (based on simulations run)
nsteps = 2000
ntraj = 24
time_step = 0.5 # ns

# statistics
nboot = 200
trials = 200

msd = []
for i in range(200):

	walker = CTRW(nsteps, ntraj, hop_dist='fbm', dwell_dist='power', hop_sigma=sigma, alpha=alpha, dt=time_step, H=hurst)

	walker.generate_trajectories(fixed_time=True)
	walker.calculate_msd(ensemble=True)
	msd.append(walker.msd.mean(axis=0)[-1])

CI = stats.confidence_interval(msd, 95)
print(np.mean(msd) - CI[1], CI[0] + np.mean(msd))
#plt.hist(msd)
#plt.show()	

walker = CTRW(nsteps, ntraj, hop_dist='fbm', dwell_dist='power', hop_sigma=sigma, alpha=alpha, dt=time_step, H=hurst)

walker.generate_trajectories(fixed_time=True)
walker.calculate_msd(ensemble=True)

walker.bootstrap_msd(nboot=nboot, fit_power_law=True)
walker.plot_msd(plot_power_law=True, show=True)


