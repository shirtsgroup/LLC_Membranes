#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import names as na

names = np.array(["ACH", "ACN", "ATO", "BUT", "DMF", "DMP", "DMS", "EAC", "ETH", "GCL", "GLY", "MET", "PCB", "PG", "PR", "RIB", "SOH", "TET", "THF", "URE"], dtype=object)

msd = np.array([66.49, 53.68, 51.85, 54.33, 54.47, 50.14, 41.83, 48.27, 51.09, 64.65, 35.58, 67.33, 52.29, 58.27, 61.91, 36.15, 46.90, 51.75, 52.86, 54.71])

error_lower = -np.array([64.40, 51.66, 50.30, 51.00, 51.95, 48.12, 40.66, 46.66, 49.00, 61.86, 34.41, 65.00, 50.18, 56.11, 59.44, 34.42, 44.90, 49.82, 50.98, 52.80]) + msd

error_upper = np.array([68.66, 56.04, 53.43, 55.94, 56.59, 52.19, 43.24, 49.81, 52.89, 66.71, 36.90, 70.00, 54.45, 60.28, 63.85, 38.01, 48.61, 53.53, 54.69, 56.70]) - msd

ordered = np.argsort(msd)[::-1]
labels = np.array([na.res_to_name[r] for r in names])[ordered]

fig, ax = plt.subplots(figsize=(12, 7))
bar_width = 0.75
index = np.arange(len(names))
ax.tick_params(labelsize=14)
ax.set_xticks(index)
ax.set_xticklabels(labels, fontsize=14)
plt.xticks(rotation=90)
plt.bar(index, msd[ordered], yerr=(error_lower[ordered], error_upper[ordered]))
plt.tight_layout()
plt.savefig('water_msd.pdf')
plt.show()
