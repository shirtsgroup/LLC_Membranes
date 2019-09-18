#!/usr/bin/env python

""" Demonstrate that the autocorrelation function for FLM has the same structure
as FBM by showing that the autocovariance scales as C * A(tau) where A is the
autocorrelation function.
"""

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.timeseries.fractional_levy_motion import FLM
import fbm
from LLC_Membranes.llclib import timeseries

H = 0.35
alpha = 1.4
nrealizations = 100 
low = 8
high = 12
C = np.zeros([high - low])


for i, j in enumerate(range(low, high)):

	l = 2**j
	flm = FLM(H, alpha, N=l, M=6000, scale=1)
	flm.generate_realizations(nrealizations)
	flm.autocovariance()
	C[i] = flm.autocov[:, 0].mean()
	flm.plot_autocorrelation(show=False, overlay=True, label='N = %d' % l, max_k=10)

acf_fbm = np.zeros([nrealizations, 9998])
for i in range(nrealizations):

	FBM = fbm.FBM(10000, H).fgn()
	acf_fbm[i, :] = timeseries.acf(FBM)

timeseries.plot_autocorrelation(acf_fbm, bootstrap=True, max_k=10, label='FBM', color='black', overlay=True)
#plt.plot(acf_fbm[:, :10].mean(), '--', lw=2, label='FBM', color='black')
plt.ylim(-0.5, 0.5)

plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('flm_autocovariance.pdf')

fig, ax = plt.subplots()
ax.plot(2 **  np.arange(low, high), C, 'o')
ax.plot(2 ** np.arange(low, high), C)
ax.set_xscale('log', basex=2)
ax.set_xlabel('Number of realization per trajectory', fontsize=14)
ax.set_ylabel('Normalization constant', fontsize=14)
plt.tight_layout()

plt.show()
