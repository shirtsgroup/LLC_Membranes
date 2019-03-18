#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw, topology
from LLC_Membranes.analysis.hbonds import System, Residue
import names
from matplotlib.patches import Patch 
import itertools

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
solutes = ["GCL", "PG", "GLY", "TET", "RIB"]
resnames = [names.res_to_name[i] for i in solutes]

d = [.3, .35, .4]
angles = [25, 30, 35]
combos = list(itertools.product(d, angles))  # all combinations of d and angles
print(combos)
for i, n in enumerate(resnames):
    if len(n.split()) == 2:
        resnames[i] = n.split()[0] + '\n' + n.split()[1]

nsolutes = 24

n = []
n_std = []
start = 0

fig, ax = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(14, 10))
bar_width = 0.2
opacity = 0.8
loc = 0

xticks = []
nhbonds = np.zeros([nsolutes], dtype=int)

for j in combos:
    yn = d.index(j[0])
    xn = angles.index(j[1])
    print(j)
    loc = 0
    for i in solutes:

        single = 0
        double = 0
        triple = 0
        quadruple = 0

        #sys = file_rw.load_object('%s/%s/10wt/hbonds.pl' % (path, i))
        sys = System('%s/%s/10wt/PR_nojump.xtc' % (path, i), '%s/%s/10wt/berendsen.gro' % (path, i), '%s/%s/10wt/topol.top' % (path, i))
        sys.set_eligible(i, 'all')
        sys.set_eligible('HII', ['O', 'O1', 'O2', 'O3', 'O4'])
        sys.identify_hbonds(j[0], j[1])
        file_rw.save_object(sys, '%s/%s/10wt/hbonds_%s_%s.pl' %(path, i, j[0], j[1]))        

        res_numbers = sys.number_residues(i)[0]
        for a in sys.hbonds:
            numbers = [res_numbers[x] for x in a[0]]
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

        start_loc = loc
        if single != 0:
            ax[xn, yn].bar(loc, single, bar_width, color='xkcd:blue', alpha=opacity)
        loc += bar_width
        if double != 0:
            ax[xn, yn].bar(loc, double, bar_width, color='xkcd:orange', alpha=opacity)
        loc += bar_width
            #print(single / double)
        if triple > 0.1:
            ax[xn, yn].bar(loc, triple, bar_width, color='xkcd:green', alpha=opacity)
        loc += bar_width
            #print(double / triple)
        if quadruple > 0.1:
            ax[xn, yn].bar(loc, quadruple, bar_width, color='xkcd:red', alpha=opacity)
        loc += bar_width
            #print(triple / quadruple)

        xticks.append((-bar_width / 2) + start_loc + (loc - start_loc) / 2)
        loc += bar_width

custom = [Patch(facecolor='xkcd:blue', alpha=opacity),
          Patch(facecolor='xkcd:orange', alpha=opacity),
          Patch(facecolor='xkcd:green', alpha=opacity),
          Patch(facecolor='xkcd:red', alpha=opacity)]

#ax[0, 2].legend(custom, ['n = 1', 'n = 2', 'n = 3', 'n = 4'], fontsize=10)

ax[0, 0].set_ylabel('$\Theta$ = 25$\degree$', fontsize=14)
ax[1, 0].set_ylabel('$\Theta$ = 30$\degree$', fontsize=14)
ax[2, 0].set_ylabel('$\Theta$ = 35$\degree$', fontsize=14)

ax[2, 0].set_xlabel('d = 3.0 $\AA$', fontsize=14)
ax[2, 1].set_xlabel('d = 3.5 $\AA$', fontsize=14)
ax[2, 2].set_xlabel('d = 4.0 $\AA$', fontsize=14)

ax[1, 1].set_xticks(xticks)
ax[1, 1].set_xticklabels(resnames, fontsize=10)
ax[1, 1].tick_params(labelsize=10)
#plt.xticks(rotation=)
#plt.ylabel('Hydrogen bond interactions per frame', fontsize=14)
plt.tight_layout()
plt.savefig('hbond_sensitivity.pdf')
plt.show()
