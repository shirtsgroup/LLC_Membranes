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

    single = 0
    double = 0
    triple = 0
    quadruple = 0

    sys = file_rw.load_object('%s/%s/10wt/hbonds.pl' % (path, i))
    res_numbers = sys.number_residues(i)[0]
    for a in sys.hbonds:
        numbers = [res_numbers[j] for j in a[0]]
        for k in numbers:
            nhbonds[k] += 1
        single += np.count_nonzero(nhbonds == 1)
        double += np.count_nonzero(nhbonds == 2)
        triple += np.count_nonzero(nhbonds == 3)
        quadruple += np.count_nonzero(nhbonds == 4)

        nhbonds -= nhbonds

    single /= sys.t.n_frames
    double /= sys.t.n_frames
    triple /= sys.t.n_frames
    quadruple /= sys.t.n_frames
    total = single + 2*double + 3*triple + 4*quadruple

    start_loc = loc

    ax.bar(loc, total, bar_width, color='white', hatch='/', alpha=opacity, edgecolor='black', label='Total')
    loc += bar_width

    if single != 0:
        ax2.bar(loc, single, bar_width, color='xkcd:blue', alpha=opacity)
    loc += bar_width
    if double != 0:
        ax2.bar(loc, double, bar_width, color='xkcd:orange', alpha=opacity)
    loc += bar_width
        #print(single / double)
    if triple > 0.1:
        ax2.bar(loc, triple, bar_width, color='xkcd:green', alpha=opacity)
    loc += bar_width
        #print(double / triple)
    if quadruple > 0.1:
        ax2.bar(loc, quadruple, bar_width, color='xkcd:red', alpha=opacity)
    loc += bar_width
        #print(triple / quadruple)

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
ax2.set_ylabel('Simultaneious Hydrogen Bond \n Interactions Per Frame', fontsize=14)
ax.set_ylabel('Total Hydrogen Bond \n Interactions Per Frame', fontsize=14)
ax2.tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('multi_hbonds.pdf')
plt.show()
