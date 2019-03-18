#!/usr/bin/env python

from LLC_Membranes.analysis.spline import Spline
from LLC_Membranes.llclib import file_rw
import numpy as np

wt = 5 
res = ['ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE']

path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11'

tortuosity = []
for i, r in enumerate(res):
	t = file_rw.load_object('%s/%s/%swt/tortuosity.pl' % (path, r, wt))
	tortuosity += list(t.flatten())

print('Average Tortuosity: %.2f +/- %.2f' % (np.mean(tortuosity), np.std(tortuosity)))
