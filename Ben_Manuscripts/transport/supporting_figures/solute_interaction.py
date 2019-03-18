#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import names

resnames = ['ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE']

labels = [names.res_to_name[i] for i in resnames]

avg_coord = [0.00120833333, 0.00091666666666, .000008333333, 0, .0000416666, 0, 0, 0.000375, 0.0004583333, 0.00054166666, 0, 0.00810441666, .0000416666666666, 0, 0, 0, 0.000125, 0, 0, 0.0020416666] 

index = np.arange(len(labels))

fig, ax = plt.subplots(figsize=(6, 7))

bar_width = 0.5
opacity = 0.7
ax.bar(index, avg_coord, bar_width, alpha=opacity)

ax.set_ylabel('Average number of\ncoordinated same-solutes', fontsize=14)
ax.tick_params(labelsize=14)
ax.set_xticks(index)
ax.set_xticklabels(labels, fontsize=14)
plt.xticks(rotation=90)

plt.tight_layout()
plt.savefig('solute_interaction.pdf')
plt.show()

