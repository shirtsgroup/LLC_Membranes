#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import names

resnames = ['MET', 'ETH', 'PR', 'BUT']

# incomplete trajectories and stricter geometric criteria
#ratio = [2.3, 3.1, 3.47, 5.51]
#total_hbonds = [8819, 9762, 8857, 8651]
ratio = [1.21781242618, 1.62801302932, 1.98014411529, 2.65940967134]
frames = 2000
#total_hbonds = 100 * np.array([18776, 20170, 18661, 17481]) / (frames * 24)
total_hbonds = np.array([43.2, 42.2, 39.5, 38.9])

index = np.arange(len(ratio))
opacity=0.75
bar_width = 0.4

fig, ax = plt.subplots()
 
ax.bar(index, ratio, bar_width, color='xkcd:blue', alpha=opacity)
ax.tick_params(axis='both', labelsize=14)
ax.set_ylabel('Carboxylate : Ether Hydrogen Bonds', fontsize=14, color='xkcd:blue')
ax.tick_params(axis='y', labelcolor='xkcd:blue')
ax2 = ax.twinx()
ax2.bar(index + bar_width, total_hbonds, bar_width, color='xkcd:red', alpha=opacity)
ax2.tick_params(axis='both', labelsize=14, labelcolor='xkcd:red')
ax2.set_ylabel('% Solutes Hydrogen Bonded Per Frame', fontsize=14, color='xkcd:red')

ax.set_xticks(index + bar_width/2)
ax.set_xticklabels([names.res_to_name[r] for r in resnames], fontsize=14)
plt.tight_layout()
plt.savefig('simple_alcohol_hbonds.pdf')
plt.show()
