#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

URE = {25: [99580.12, 11229.48], 30: [165503.50, 18697.03], 35: [259302.96, 50973.57]}

x = list(URE.keys())
y = np.array([i for i in URE.values()])

plt.errorbar(x, y[:, 0] / 1000, yerr=y[:, 1] / 1000, marker="o", lw=2, elinewidth=2, capsize=10, capthick=2)
plt.tick_params(labelsize=14)
plt.xlabel('Pore Length (nm)', fontsize=14)
plt.ylabel('Mean First Passage Time ($\mu$s)', fontsize=14)
plt.tight_layout()
plt.show()
