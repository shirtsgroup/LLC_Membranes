#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

L = [1, 5, 10, 25, 50, 100]
#rpi = [74.02, __, 72.76, __, 76.16]

v = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1])
rpi_v = np.array([100, 72.76, 56.99, 42.58, 31.65, 25.26, 13.59, 8])

correlation = np.load('correlation.npz')
no_correlation = np.load('nocorrelation.npz')

#plt.figure()
#plt.plot(v, rpi_v)
#plt.xlabel('Variance in scatterer position')
#plt.ylabel('Intensity')
#plt.tight_layout()
#plt.savefig('correlation_decay.png')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(no_correlation['freq_z'], no_correlation['slice'], label='No Correlation', linewidth=2)
ax.plot(correlation['freq_z'], correlation['slice'], label='Correlation', linewidth=2)
left, bottom, width, height = [0.38, 0.4, 0.35, 0.35]
ax2 = fig.add_axes([left, bottom, width, height])
ax2.plot(no_correlation['freq_z'], no_correlation['slice'], label='No Correlation', linewidth=2)
ax2.plot(correlation['freq_z'], correlation['slice'], label='Correlated', linewidth=2)
ax2.set_ylim(0, 20)
ax2.set_xlim(1.2, 1.6)
#plt.ylim(0, 30)
#plt.xlim(0.5, 4)
ax.legend(loc=9, fontsize=14)
ax.set_xlabel('$q_z\ (\AA^{-1})$', fontsize=14)
ax.set_ylabel('Intensity', fontsize=14)
ax.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
#plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('../../structure_paper/figures/sf_qz_correlation.png')
plt.show()
