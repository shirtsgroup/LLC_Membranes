#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
from LLC_Membranes.analysis.ztrace import ZTrace
from LLC_Membranes.llclib import file_rw

solute = "URE"
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt" % solute
traj = "5ms_nojump.xtc"
gro = "em.gro"
load = True 

ndx = [4, 7, 18]

if load: 
	trace = file_rw.load_object('%s_trace.pl' % solute)
else:
	trace = ZTrace("%s/%s" % (path, traj), "%s/%s" %(path, gro), solute, 'HII', 'z')
	trace.locate_pore_centers(save=True, savename="%s/spline.pl" % path)
	trace.radial_distances()
	trace.t = None
	file_rw.save_object(trace, '%s_trace.pl' % solute)

trace.plot_trace(ndx, cmax=1.5, savename='%s_trajectories.pdf' % solute)
#for i in range(24):
#	trace.plot_trace(i, cmax=1.5)
#t = md.load("%s/%s" % (path, traj), top="%s/%s" % (path, gro))
#keep = [a.index for a in t.topology.atoms if a.residue.name == 
