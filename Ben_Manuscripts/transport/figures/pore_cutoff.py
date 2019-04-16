#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

cut = [0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80]

diff = [44.59, 45.17, 45.47, 45.25, 44.09, 43.12, 44.98, 46.30, 52.64, 53.22, 53.86, 55.23, 58.12, 59.39, 59.26, 58.71, 58.05, 57.68, 58.15]

plt.plot(cut, diff, linewidth=2)
#plt.plot([0.73, 0.73], [min(diff), max(diff)], '--', color='black', linewidth=2)
plt.tick_params(labelsize=14)
plt.xlabel('Pore region cut-off (nm)', fontsize=14)
plt.ylabel('Hop length in pore relative to tails (%)', fontsize=14)
plt.tight_layout()
plt.savefig('pore_cutoff.pdf')
plt.show()

