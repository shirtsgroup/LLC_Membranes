#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import names

names = [names.res_to_name[i] for i in ['ACH', 'URE', 'ACN', 'ATO']]

hg = (100 / 24) * np.array([12.7, 1.84, 1.12, 0])
donors = (100 / 24) * np.array([6.6, 2.6, 1.2, 0])
acc = (100 / 24) * np.array([1.45, 2.2, 1.6, 0.5])
na = (100 / 24) * np.array([.28, .52, .45, .26]) * 24

index = np.arange(len(names))

fig, ax = plt.subplots()

bar_width = 0.2

ax.set_ylabel('% Solutes Involved Per Frame', fontsize=14)
ax.set_xlabel('   ')
ax.tick_params(labelsize=14)
ax.set_xticks(index + 2*bar_width)
ax.set_xticklabels(names, fontsize=14)

opacity = 0.7

ax.bar(index, hg, bar_width, alpha=opacity, label='hbonds w/ head groups')
ax.bar(index + bar_width, donors, bar_width, alpha=opacity, label='hbonds donated to water')
ax.bar(index + 2*bar_width, acc, bar_width, alpha=opacity, label='hbonds accepted from water')
ax.bar(index + 3*bar_width, na, bar_width, alpha=opacity, label='coordinated to sodium')

plt.ylim(0, 80)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('ketone_hbonds.pdf')
plt.show()
