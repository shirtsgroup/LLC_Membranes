#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw, topology
from LLC_Membranes.analysis.hbonds import System, Residue
import names
from matplotlib.patches import Patch 

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
solutes = ["GCL", "PG", "GLY", "TET", "RIB"]
resnames = [names.res_to_name[i] for i in solutes]
for i, n in enumerate(resnames):
    if len(n.split()) == 2:
        resnames[i] = n.split()[0] + '\n' + n.split()[1]
nsolutes = 24
nboot = 200  # number of bootstrap trials

n = []
n_std = []
start = 0

fig, ax = plt.subplots(figsize=(7, 5))
ax2 = ax.twinx()
bar_width = 0.2
opacity = 0.8
loc = 0

xticks = []
nhbonds = np.zeros([nsolutes], dtype=int)
for i in solutes:

    sys = file_rw.load_object('%s/%s/10wt/hbonds.pl' % (path, i))
    res_numbers = sys.number_residues(i)[0]

    single = np.zeros([sys.t.n_frames, nsolutes], dtype=bool)
    double = np.zeros_like(single) 
    triple = np.zeros_like(single)
    quadruple = np.zeros_like(single)

    for t, a in enumerate(sys.hbonds):
        nhbonds = np.zeros([nsolutes], dtype=int)
        numbers = [res_numbers[j] for j in a[0]]
        for k in numbers:
            nhbonds[k] += 1

        single[t, np.where(nhbonds == 1)[0]] = True
        double[t, np.where(nhbonds == 2)[0]] = True
        triple[t, np.where(nhbonds == 3)[0]] = True
        quadruple[t, np.where(nhbonds == 4)[0]] = True

    boot = np.zeros([5, nboot])

    for b in range(nboot):
        ndx = np.random.randint(0, nsolutes, size=nsolutes)
        boot[0, b] = single[:, ndx].sum(axis=1).mean()
        boot[1, b] = double[:, ndx].sum(axis=1).mean()
        boot[2, b] = triple[:, ndx].sum(axis=1).mean()
        boot[3, b] = quadruple[:, ndx].sum(axis=1).mean()
        boot[4, b] = boot[0, b] + 2 * boot[1, b] + 3 * boot[2, b] + 4 * boot[3, b]

    total = boot[0, :].mean() + 2*boot[1, :].mean() + 3*boot[2, :].mean() + 4*boot[3, :].mean()
    start_loc = loc

    ax.bar(loc, boot[4, :].mean(), bar_width, color='white', hatch='/', alpha=opacity, edgecolor='black', label='Total', yerr=boot[4, :].std())
    loc += bar_width

    ax2.bar(loc, boot[0, :].mean(), bar_width, color='xkcd:blue', alpha=opacity, yerr=boot[0, :].std())
    loc += bar_width
    
    ax2.bar(loc, boot[1, :].mean(), bar_width, color='xkcd:orange', alpha=opacity, yerr=boot[1, :].std())
    loc += bar_width

    ax2.bar(loc, boot[2, :].mean(), bar_width, color='xkcd:green', alpha=opacity, yerr=boot[2, :].std())
    loc += bar_width

    ax2.bar(loc, boot[3, :].mean(), bar_width, color='xkcd:red', alpha=opacity, yerr=boot[3, :].std())
    loc += bar_width

    xticks.append((-bar_width / 2) + start_loc + (loc - start_loc) / 2)
    loc += bar_width

custom = [Patch(facecolor='xkcd:blue', alpha=opacity),
          Patch(facecolor='xkcd:orange', alpha=opacity),
          Patch(facecolor='xkcd:green', alpha=opacity),
          Patch(facecolor='xkcd:red', alpha=opacity),
          Patch(facecolor='white', alpha=opacity, hatch='/', edgecolor='black')]

ax2.legend(custom, ['n = 1', 'n = 2', 'n = 3', 'n = 4', 'Total Hbonds'], fontsize=14, ncol=3, bbox_to_anchor=(1, 1.25))

ax.set_xticks(xticks)
ax.set_xticklabels(resnames, fontsize=14)
ax.tick_params(labelsize=14)
#plt.xticks(rotation=)
ax2.set_ylabel('Simultaneous Hydrogen Bond \n Interactions Per Frame', fontsize=14)
ax.set_ylabel('Total Hydrogen Bond \n Interactions Per Frame', fontsize=14)
ax2.tick_params(labelsize=14)
ax2.set_ylim(0, 11.5)
plt.tight_layout()
plt.savefig('multi_hbonds.pdf')
plt.show()
