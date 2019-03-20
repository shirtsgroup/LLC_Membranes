#!/usr/bin/env python

from LLC_Membranes.analysis.coordination_number import System
from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt

head_groups = True
wts = [5, 10] 
#res = ['ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE']

res = 'MET'
path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11'

if head_groups:
	pickle = 'hg_coordination.pl'
	savename = 'Na_hg_coordination.pdf'
	coordinated = 'Carboxylate Groups'
else:
	pickle = 'coordination.pl'
	savename = 'Na_water_coordination.pdf'
	coordinated = 'Water Molecules'

for wt in wts:

	c = file_rw.load_object('%s/%s/%dwt/%s' % (path, res, wt, pickle))

	n = [(c.distances[t] !=0).sum(1).mean() for t in range(len(c.distances))]

	plt.plot(c.t.time / 1000, n, label='%d wt %% water' % wt, linewidth=2)

plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.xlabel('Time (ns)', fontsize=14)
plt.ylabel('Average %s \n Coordinated to Sodium' % coordinated, fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
#plt.savefig(savename)
plt.show()

