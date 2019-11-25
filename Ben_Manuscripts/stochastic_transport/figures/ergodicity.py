#!/usr/bin/env python

from LLC_Membranes.analysis.msd import Diffusivity
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import numpy as np

res = "URE"
gro = "em.gro"
traj = "5ms_nojump.xtc"
directory = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt" % res
first_frame = {"URE": 2000, "GCL": 2400}
nboot = 200
fracshow=0.4
dt = 0.5
load = True
stationarity = False 

if load:
	msd = file_rw.load_object('msd.pl')
else:
	msd = Diffusivity('%s/%s' % (directory,traj), '%s/%s' % (directory, gro), 'z', begin=first_frame[res], residue=res)
	file_rw.save_object(msd, 'msd.pl')


if stationarity:
	p = msd.stationarity(axis=2)
	print('Probability that time series is non-stationary:')
	for i, pvalue in enumerate(p):
		print('Trajectory %d: %.2f' % (i + 1, pvalue))

msd.calculate()
#eb = msd.ergodicity_breaking_parameter()
#plt.plot(eb)
#plt.show()
#exit()

msd.bootstrap(nboot, fit_line=False)

end = int(fracshow * msd.time.size)
t = np.arange(end) * dt

plt.plot(t, msd.MSD_average[:end], label='Time Averaged MSD')
plt.fill_between(t, msd.MSD_average[:end] + msd.limits[0, :end], msd.MSD_average[:end] - msd.limits[1, :end], alpha=0.7)

msd.calculate(ensemble=True)
msd.bootstrap(nboot, fit_line=False)
plt.plot(t, msd.MSD_average[:end], label='Ensemble Averaged MSD')
plt.fill_between(t, msd.MSD_average[:end] + msd.limits[0, :end], msd.MSD_average[:end] - msd.limits[1, :end], alpha=0.7)

plt.legend(fontsize=14)
plt.show()
#reduce object size and save
#MD_MSD.t = None
#MD_MSD.com = None
#file_rw.save_object(MD_MSD, '%s_msd.pl' % res)


