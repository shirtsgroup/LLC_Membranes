#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

pd = [0.597, .668, 0.843, 0.926, 1.092]
pd_error = [0.008, .003, 0.02, 0.015, 0.007]

s = [0.578, .666, 0.813, 0.846, 0.969]
s_error = [0.009, .007, 0.013, 0.012, 0.023]

s_disordered = [.632, .717]
s_disordered_error = [.004, 0.01]

pd_disordered = [.638, .706]
pd_disordered_error = [.011, .013]

width = 0.2

pos = np.arange(5)
pos = np.array([p/2 for p in pos])
pos[2:] += 2*width
pos[3:] += 2*width

fig, ax = plt.subplots(figsize=(10, 5))

ax.bar(pos, pd, width, yerr=pd_error, label='Sandwiched (d=3.7)')
ax.bar([p + width for p in pos], s, width, yerr=s_error, label='Parallel Displaced (d=3.7)')
ax.bar([p + 2*width for p in pos][1:3], s_disordered, width, label='Sandwiched (d=5)', yerr=s_disordered_error)
ax.bar([p + 3*width for p in pos][1:3], pd_disordered, width, label='Parallel Displaced (d=5)', yerr=pd_disordered_error)

tick_locations = [p + .5*width for p in pos]
tick_locations[1] += width
tick_locations[2] += width
ax.set_xticks(tick_locations)
ax.set_xticklabels([4, 5, 6, 7, 8], fontsize=14)
ax.set_ylabel('Pore Radius (nm)', fontsize=14)
ax.set_xlabel('Number of columns per pore', fontsize=14)
plt.tick_params(axis='both', labelsize=14)
#ax.plot([-0.20, 3.2], [4.12, 4.12], "k--", label='Experiment')
#plt.legend(['Sandwiched', 'Parallel Displaced', 'Experiment'],fontsize=14)
plt.legend(fontsize=14)
plt.ylim(0.5, 1.2)
plt.xlim(-0.2, 3.2)
plt.tight_layout()
plt.savefig('pore_radius.pdf')
plt.show()

plt.show()
