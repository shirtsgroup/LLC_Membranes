#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw, topology
from LLC_Membranes.analysis.hbonds import System, Topology, Residue
import names
from matplotlib.patches import Patch 

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"

regular = ["GCL", "GLY"]
analog = ["SOH", "DMP"]

#resnames = [names.res_to_name[i] for i in solutes]
#for i, n in enumerate(resnames):
#    if len(n.split()) == 2:
#        resnames[i] = n.split()[0] + '\n' + n.split()[1]

nsolutes = 24

n = []
n_std = []
start = 0

fig, ax = plt.subplots(figsize=(10, 5))
bar_width = 0.2
opacity = 0.8
loc = 0

xticks = []
nhbonds = np.zeros([nsolutes], dtype=int)
for ndx, i in enumerate(regular):

    print(analog[ndx], i)

    reg = file_rw.load_object('%s/%s/10wt/hbonds.pl' % (path, i))
    ana = file_rw.load_object('%s/%s/10wt/hbonds.pl' %(path, analog[ndx]))
    
    total_reg = 0
    for k in reg.hbonds:
        total_reg += k.shape[1]

    total_analog = 0
    for k in ana.hbonds:
        total_analog += k.shape[1]

    print(total_reg / total_analog)

#custom = [Patch(facecolor='xkcd:blue', alpha=opacity),
#          Patch(facecolor='xkcd:orange', alpha=opacity),
#          Patch(facecolor='xkcd:green', alpha=opacity),
#          Patch(facecolor='xkcd:red', alpha=opacity)]

#ax.legend(custom, ['n = 1', 'n = 2', 'n = 3', 'n = 4'], fontsize=14)

ax.set_xticks(xticks)
ax.set_xticklabels(resnames, fontsize=14)
ax.tick_params(labelsize=14)
#plt.xticks(rotation=)
ax.set_ylabel('Hydrogen bond interactions', fontsize=14)
plt.tight_layout()
plt.savefig('multi_hbonds.pdf')
plt.show()
