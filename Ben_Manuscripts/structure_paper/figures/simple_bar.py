#!/usr/bin/env python
import matplotlib.pyplot as plt

# dbwl bar plot
x = [0.375, 0.625]
y = [4.43, 4.39]

width = 0.25

fig, ax = plt.subplots(figsize=(10,5))

plt.bar(x[0], y[0], width=width, label='xlabel')
plt.bar(x[1], y[1], width=width, label='label2')

ax.set_xticks(x)
ax.set_xticklabels(['Crosslinked', 'Uncrosslinked'],fontsize=14)

ax.set_aspect('2')
plt.tick_params(axis='both', labelsize=14)

plt.ylim(4.25, 4.5)
plt.ylabel('Distance between layers ($\AA$)', fontsize=14)
plt.tight_layout()
plt.savefig('dbwl_xlink.png')

plt.show()
exit()
# p2p bar plot

x = [0.375, 0.625]
y = [4.19, 4.23]
y_err = [0.04, 0.04]

width = 0.25

fig, ax = plt.subplots(figsize=(10,5))

plt.bar(x[0], y[0], width=width, yerr=y_err[0])
plt.bar(x[1], y[1], width=width, yerr=y_err[1])

ax.set_xticks(x)
ax.set_xticklabels(['Crosslinked', 'Uncrosslinked'],fontsize=14)

ax.set_aspect('1.666666')
plt.tick_params(axis='both', labelsize=14)

plt.ylim(4, 4.3)
plt.ylabel('Pore-to-pore distance (nm)', fontsize=14)
plt.tight_layout()
plt.savefig('p2p_xlink.png')

plt.show()
