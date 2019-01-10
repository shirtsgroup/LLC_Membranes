#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

"""
Hurst parameter estimation
"""

def autocovariance(H, length):

	k = np.arange(length)

	return 0.5*(np.abs(k - 1)**(2*H) - 2*k**(2*H) + (k + 1)**(2*H))

L = 10
h = np.linspace(0.05, 0.5, 100)

simulated_h = np.linspace(0.05, 0.5, 10)
simulated_autocorrelation = [-.464, -.424, -.384, -.339, -.292, -.243, -.188, -.127, -.067, -.0014]  # fbmsim.py -H $h -acf -nftraj 1000
simulated_autocov = [-.4643, -.4252, -.3842, -.3403 ,-.2935, -.2404 ,-.1868, -.1286, -0.06729, 0.00024]  # fbmsim.py -H $h -acov -nftraj 1000

alpha = 1

plt.plot(h, [autocovariance(x, L)[1] for x in h], linewidth=2, alpha=alpha, label='Analytical')
plt.plot(simulated_h, simulated_autocorrelation, '--', linewidth=2, alpha=alpha, label='Simulated Autocorrelation')
plt.plot(simulated_h, simulated_autocov, 'o', linewidth=2, alpha=alpha, label='Simulated Autocovariance')

plt.xlabel('Hurst parameter', fontsize=14)
plt.ylabel('Autocovariance / Autocorrelation \n of first timestep (k = 1)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.legend()

plt.show()
