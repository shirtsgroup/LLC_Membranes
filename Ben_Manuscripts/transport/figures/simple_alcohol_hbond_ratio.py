#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import names

resnames = ['MET', 'ETH', 'PR', 'BUT']

ratio = [2.3, 3.1, 3.47, 5.51]
total_hbonds = [8819, 9762, 8857, 8651]

index = np.arange(len(ratio))
opacity=0.9
bar_width = 0.4

fig, ax = plt.subplots()
 
ax.bar(index, ratio, bar_width, color='xkcd:blue', alpha=opacity)
ax.tick_params(axis='both', labelsize=14)
ax.set_ylabel('Carboxylate : Ether Hydrogen Bonds', fontsize=14, color='xkcd:blue')
ax.tick_params(axis='y', labelcolor='xkcd:blue')
ax2 = ax.twinx()
ax2.bar(index + bar_width, total_hbonds, bar_width, color='xkcd:red', alpha=opacity)
ax2.tick_params(axis='both', labelsize=14, labelcolor='xkcd:red')
ax2.set_ylabel('Total Hydrogen Bond Occurences', fontsize=14, color='xkcd:red')

ax.set_xticks(index + bar_width/2)
ax.set_xticklabels([names.res_to_name[r] for r in resnames], fontsize=14)
plt.tight_layout()
plt.savefig('simple_alcohol_hbonds.pdf')
plt.show()
