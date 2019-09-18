#!/usr/bin/env python

from LLC_Membranes.timeseries.fractional_levy_motion import FLM
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import levy_stable

H = 0.35
alpha = 1.4
l = 2**12
nrealizations = 1
opacity=0.7
truncate = 1
bins = 100 
plot_range=(-0.6, 0.6)

flm = FLM(H, alpha, M=2, N=l, scale=0.035, truncate=0.5, correct_truncation=False)
flm.generate_realizations(nrealizations)
plt.hist(flm.noise.flatten(), bins=bins, range=plot_range, alpha=opacity, density=True, label='Uncorrected Truncation Parameter\nMax Magnitude: %.2f' % np.abs(flm.noise.flatten()).max())

flm = FLM(H, alpha, M=2, N=l, scale=0.035, truncate=0.5, correct_truncation=True)
flm.generate_realizations(nrealizations)
plt.hist(flm.noise.flatten(), bins=bins, range=plot_range, alpha=opacity, density=True, label='Corrected Truncation Parameter\nMax Magnitude: %.2f' % np.abs(flm.noise).max())

x = np.linspace(plot_range[0], plot_range[1], 100)
plt.plot(x, levy_stable.pdf(x, alpha=alpha, loc=0, scale=0.035, beta=0), '--', lw=2, color='black')

plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=12)
plt.ylabel('Probability Density', fontsize=14)
plt.tight_layout()
plt.savefig('truncation_correction.pdf')
plt.show()
