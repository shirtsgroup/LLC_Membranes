#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import timeseries

ntraj = 10 
npoints = 10000
fracshow = 0.2
low_T = 100

# make Brownian trajectories

x = np.zeros([npoints, ntraj])
for t in range(ntraj):
	x[:, t] = np.cumsum(np.random.normal(size=npoints))

msd = timeseries.msd(x[..., np.newaxis], 0)
limits = timeseries.bootstrap_msd(msd, 200)

msd_mean = msd.mean(axis=1)
msd_mean_squared = np.square(msd_mean)

mean_msd_variance = np.zeros_like(msd_mean)
for i in range(low_T, npoints):
	mean_msd_variance[i] = np.square(msd[i, :]).mean()

eb = (mean_msd_variance - msd_mean_squared) / msd_mean_squared

t = np.arange(npoints)
end = int(fracshow * npoints)

plt.plot(t[:end], msd_mean[:end])
plt.fill_between(t[:end], msd_mean[:end] + limits[0, :end], msd_mean[:end] - limits[1, :end], alpha=0.7)

plt.figure()
plt.plot(t[low_T:], eb[low_T:])
plt.show()

