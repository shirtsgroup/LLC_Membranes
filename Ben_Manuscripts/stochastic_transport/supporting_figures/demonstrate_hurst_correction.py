#!/usr/bin/env python

from LLC_Membranes.timeseries.fractional_levy_motion import FLM
import numpy as np
import matplotlib.pyplot as plt

H = 0.35
alpha = 1.4
l = 2**12
nrealizations = 100
max_k = 10

flm = FLM(H, alpha, M=2, N=l, scale=0.035, correct_hurst=False)
flm.generate_realizations(nrealizations)
flm.autocorrelation()
h_estimate = np.log(2 * flm.acf[:, 1].mean() + 2) / (2 * np.log(2))
flm.plot_autocorrelation(label='Uncorrected Input H \n Estimated H: %.2f' % h_estimate, max_k=max_k)

flm = FLM(H, alpha, M=2, N=l, scale=0.035, correct_hurst=True)
flm.generate_realizations(nrealizations)
flm.autocorrelation()
h_estimate = np.log(2 * flm.acf[:, 1].mean() + 2) / (2 * np.log(2))
flm.plot_autocorrelation(label='Corrected Input H \n Estimated H: %.2f' % h_estimate, overlay=True, max_k=max_k)

plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('hurst_correction.pdf')
plt.show()
