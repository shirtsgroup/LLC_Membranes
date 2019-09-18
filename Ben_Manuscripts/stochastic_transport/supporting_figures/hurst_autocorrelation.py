#!/usr/bin/env python

""" Plot the analytical autocorrelation function of a fractional process for
various value of the hurst parameter
"""

import numpy as np
from LLC_Membranes.llclib.fitting_functions import hurst_autocovariance as hurst
import matplotlib.pyplot as plt

H = np.linspace(0.1, 0.9, 9)
k = np.arange(1, 15)

plt.figure(figsize=(7, 5))

for h in H:
	plt.plot(k, hurst(k, h), lw=2, label='H = %.1f' % h)

plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.xlabel('Lag (k)', fontsize=14)
plt.ylabel('Autocorrelation', fontsize=14)
plt.legend(fontsize=14, ncol=3)
plt.tight_layout()
plt.savefig('hurst_autocorrelation.pdf')
plt.show()
