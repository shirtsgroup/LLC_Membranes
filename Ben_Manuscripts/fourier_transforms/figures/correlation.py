#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

L = [1, 5, 10, 25, 50, 100]
#rpi = [74.02, __, 72.76, __, 76.16]

v = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1])
rpi_v = np.array([100, 72.76, 56.99, 42.58, 31.65, 25.26, 13.59, 8])

correlation = np.load('correlation.npz')
no_correlation = np.load('no_correlation.npz')

plt.figure()
plt.plot(v, rpi_v)
plt.xlabel('Variance in scatterer position')
plt.ylabel('Intensity')
plt.tight_layout()
plt.savefig('correlation_decay.png')

plt.figure()
plt.plot(no_correlation['freq_z'], no_correlation['slice'], label='No Correlation')
plt.plot(correlation['freq_z'], correlation['slice'], label='Correlated')
plt.ylim(0, 30)
plt.xlim(0.5, 4)
plt.legend(loc=2)
plt.xlabel('$q_z~(\AA^{-1})$')
plt.ylabel('Intensity')
plt.tight_layout()
plt.savefig('correlation.png')
plt.show()
