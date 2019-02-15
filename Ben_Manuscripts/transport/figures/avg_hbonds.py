#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw, topology
from LLC_Membranes.analysis.hbonds import System, Topology, Residue
import names

## Former avg_hbonds.py moved to hbonds.py ##

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
#solutes = ["GLY", "GCL", "PG"]
solutes = ["TET", "RIB"]
solutes = ["GCL", "SOH", "GLY", "DMP"]
solutes = ["MET", "ETH", "PR", "BUT", "GCL", "PG", "GLY", "TET", "RIB"]
resnames = [names.res_to_name[i] for i in solutes]
nsolutes = 24
npts_ma = 8  # number of points to average over in moving average

n = []
n_std = []
start = 0 

for i in solutes:
    sys = file_rw.load_object('%s/%s/10wt/hbonds.pl' % (path, i))
    hbonds = [a.shape[1] / nsolutes for a in sys.hbonds]
    n.append(np.mean(hbonds))
    n_std.append(np.std(hbonds))

plt.bar(resnames, n) #, yerr=n_std) do error bars make sense for high fluctuations?
plt.xticks(rotation=75)
plt.ylabel('Average hydrogen bonds \n per solute', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.show()
