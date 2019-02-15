#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw, topology
from LLC_Membranes.analysis.hbonds import System, Topology, Residue
import names

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

fig, ax = plt.subplots(figsize=(10, 5))
bar_width = 0.2
opacity = 0.8
loc = 0

xticks = []
for nd, i in enumerate(solutes):
    n = 0

    sys = file_rw.load_object('%s/%s/10wt/hbonds.pl' % (path, i))
    res_numbers = sys.number_residues(i)[0]  # dictionary for converting index to residue number (residue 0 through residue 23)

    for a in sys.hbonds:
       for ndx, k in enumerate(a[0]):
           donor = res_numbers[k]
           try:
               acceptor = res_numbers[a[2][ndx]]
               if acceptor == donor:
                   n += 1
           except KeyError:
               pass
#        donors = np.array([res_numbers[int(j)] for j in a[0]], dtype=int)
#        acceptors = np.array([res_numbers[int(j)] for j in a[2]], dtype=int)
#        intra = (donors == acceptors)  # find cases where a donor and acceptor are on the same solute

    ax.bar(resnames[nd], n, bar_width)

#from matplotlib.patches import Patch 
#custom = [Line2D([0], [0], color='xkcd:blue', lw=4, alpha=opacity),
#          Line2D([0], [0], color='xkcd:orange', lw=4, alpha=opacity),
#          Line2D([0], [0], color='xkcd:green', lw=4, alpha=opacity),
#          Line2D([0], [0], color='xkcd:red', lw=4, alpha=opacity)]
#custom = [Patch(facecolor='xkcd:blue', alpha=opacity),
#          Patch(facecolor='xkcd:orange', alpha=opacity),
#          Patch(facecolor='xkcd:green', alpha=opacity),
#          Patch(facecolor='xkcd:red', alpha=opacity)]

#ax.legend(custom, ['n = 1', 'n = 2', 'n = 3', 'n = 4'], fontsize=14)

#ax.set_xticks(xticks)
ax.set_xticklabels(resnames, fontsize=14)
ax.tick_params(labelsize=14)
#plt.xticks(rotation=)
ax.set_ylabel('Hydrogen bond interactions', fontsize=14)
plt.tight_layout()
plt.show()
