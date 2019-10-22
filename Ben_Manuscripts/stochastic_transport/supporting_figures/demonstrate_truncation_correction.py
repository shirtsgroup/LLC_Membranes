#!/usr/bin/env python

from LLC_Membranes.timeseries.fractional_levy_motion import FLM
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import levy_stable

H = 0.37
alpha = 1.8
l = 2**12
nrealizations = 1
opacity=0.7
truncate = 1.73 
bins = 100 
plot_range=(-2, 2)
scale = 0.22

flm = FLM(H, alpha, M=2, N=l, scale=scale, truncate=truncate, correct_truncation=False)
flm.generate_realizations(nrealizations)
plt.hist(flm.noise.flatten(), bins=bins, range=plot_range, alpha=opacity, density=True, label='Uncorrected Truncation Parameter\nMax Magnitude: %.2f' % np.abs(flm.noise.flatten()).max())

flm = FLM(H, alpha, M=2, N=l, scale=scale, truncate=truncate, correct_truncation=True)
flm.generate_realizations(nrealizations)
plt.hist(flm.noise.flatten(), bins=bins, range=plot_range, alpha=opacity, density=True, label='Corrected Truncation Parameter\nMax Magnitude: %.2f' % np.abs(flm.noise).max())

x = np.linspace(plot_range[0], plot_range[1], 100)
plt.plot(x, levy_stable.pdf(x, alpha=alpha, loc=0, scale=scale, beta=0), '--', lw=2, color='black')

plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=12)
plt.ylabel('Probability Density', fontsize=14)
plt.tight_layout()
#plt.savefig('truncation_correction.pdf')
plt.show()
