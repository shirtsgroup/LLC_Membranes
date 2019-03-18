#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import names

#names = [names.res_to_name[i] for i in ['THF', 'PCB', 'EAC', 'DMF']]
names = ["Tetrahydrofuran", "Propylene\nCarbonate", "Ethyl\nAcetate", "Dimethyl\nFormamide"]
#hg = [12.7, 1.84, 1.12, 0]
#donors = [6.6, 2.6, 1.2, 0]
acc = [0.4475, 1.5, .65, .765]
na = np.array([.157, .292, .27, .28]) * 24

index = np.arange(len(names))

fig, ax = plt.subplots()

bar_width = 0.4

ax.set_ylabel('Occurences Per Frame', fontsize=14)
ax.set_xlabel('   ')
ax.tick_params(labelsize=14)
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(names, fontsize=14)

opacity = 0.7

ax.bar(index, na, bar_width, alpha=opacity, label='Coordinated to sodium')
#ax.bar(index + bar_width, donors, bar_width, alpha=opacity, label='hbonds donated to water')
ax.bar(index + bar_width, acc, bar_width, alpha=opacity, label='Hbonds accepted from water')
#ax.bar(index + 3*bar_width, na, bar_width, alpha=opacity, label='coordinated to sodium')

plt.ylim(0, 10)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('nondonor_hbonds.pdf')
plt.show()
