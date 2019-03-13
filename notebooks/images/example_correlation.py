#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

z = np.linspace(0, 4, 1000)
y_perfect = np.zeros_like(z)
y_perfect[250] = 1
y_perfect[500] = 1
y_perfect[750] = 1
y_perfect[-2] = 1

fig, ax = plt.subplots()

ax.plot(z, y_perfect, linewidth=2)
ax.set_xlabel('Distance', fontsize=14)
ax.set_ylabel('Count', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.yticks([])
plt.xticks([0, 1, 2, 3, 4])
ax.set_xticklabels(['0', 'd', '2d', '3d', '4d'])

plt.tight_layout()
plt.savefig('ordered_correlation_function.png')

fig, ax = plt.subplots()
y_disordered = np.exp(-z)*np.cos(2*np.pi*z)

mini = np.argmin(y_disordered)
shift = y_disordered[mini]
y_disordered -= shift
y_disordered[:mini] = 0


ax.plot(z, y_disordered, linewidth=2)
ax.set_ylabel('Count', fontsize=14)
ax.set_xlabel('Distance', fontsize=14)
plt.yticks([])
plt.xticks([0, 1, 2, 3, 4])
ax.set_xticklabels(['0', 'd', '2d', '3d', '4d'])
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()

plt.savefig('disordered_correlation_function.png')

plt.show()



