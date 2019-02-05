#! /usr/bin/env python

from LLC_Membranes.analysis import msd
from LLC_Membranes.llclib import file_rw
from LLC_Membranes.analysis.msd import Diffusivity
from LLC_Membranes.timeseries import forecast_ctrw
from LLC_Membranes.timeseries.forecast_ctrw import System
import matplotlib.pyplot as plt
import numpy as np

residue = 'ETH'
load = True 
path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt/' % residue

if load:
	D = file_rw.load_object('%s/msd_%s.pl' % (path, residue))
else:
	D = msd.Diffusivity('%s/PR_nojump.xtc' % path, '%s/berendsen.gro' % path, 'z', residue='%s' % residue)

	D.calculate(ensemble=False)  # do time averaged msd

	D.step_autocovariance()

	D.bootstrap(200)

	file_rw.save_object(D, '%s/msd_%s.pl' % (path, residue))

ctrw = file_rw.load_object('%s/forecast_%s.pl' % (path, residue))

small_figure_font = 21

# hop distribution
plt.figure()
hops = []
for i in ctrw.hop_lengths:
	hops += i
plt.hist(hops, bins=30, density=True)
plt.xlabel('$z$-direction Hop Length (nm)', fontsize=small_figure_font)
plt.ylabel('Frequency', fontsize=small_figure_font)
plt.gcf().get_axes()[0].tick_params(labelsize=small_figure_font)
plt.xlim(-5*np.std(hops), 5*np.std(hops))
plt.tight_layout()
plt.savefig('example_hop_dist.pdf')

# msd
plt.figure()
fracshow = 0.4
end = int(D.time.size * fracshow)
plt.plot(D.time[:end], D.MSD_average[:end])
plt.fill_between(D.time[:end], D.MSD_average[:end] + D.limits[0, :end], D.MSD_average[:end] - D.limits[1, :end], alpha=0.7)

plt.ylabel('MSD ($nm^2$)', fontsize=14)
plt.xlabel('time (ns)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('example_msd.pdf')

# ztraces
plt.figure()
np.random.seed(4)  # ethanol - 4 , methanol - 5,33, 46 (16 interesting)
trajs = np.random.randint(0, D.com.shape[1], size=3)
for i in trajs:
     plt.plot(D.time, D.com[:, i, 2], linewidth=2)

plt.ylabel('$z$ position ($nm$)', fontsize=14)
plt.xlabel('time (ns)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.tight_layout()
plt.savefig('example_ztraces.pdf')

# power law
plt.figure()
bins = plt.hist(ctrw.dwell_times, bins=25, density=True)
#lowest, highest = bins[1][0] + (bins[1][1] - bins[1][0])/2, bins[1][-1]
#xvalues = np.linspace(lowest, highest, 100)
#plt.plot(xvalues, xvalues**-(1 + np.mean(ctrw.alpha_distribution)), '--', color='black', linewidth=2)
plt.xlabel('Dwell time (ns)', fontsize=small_figure_font)
plt.ylabel('Frequency', fontsize=small_figure_font)
plt.gcf().get_axes()[0].tick_params(labelsize=small_figure_font)
plt.tight_layout()
plt.savefig('example_powerlaw.pdf')

# autocovariance
plt.figure()
plt.plot(np.arange(len(ctrw.hop_acf)), ctrw.hop_acf, linewidth=2)
plt.ylabel('Autocovariance', fontsize=small_figure_font)
plt.xlabel('k', fontsize=small_figure_font)
plt.gcf().get_axes()[0].tick_params(labelsize=small_figure_font)
plt.xlim(-1, 20)
plt.tight_layout()
plt.savefig('example_autocovariance.pdf')

plt.show()
exit()

fig, ax = plt.subplots(1, 3, figsize=(12, 5))

fracshow = 0.4
end = int(D.time.size * fracshow)
ax[0].plot(D.time[:end], D.MSD_average[:end])
ax[0].fill_between(D.time[:end], D.MSD_average[:end] + D.limits[0, :end], D.MSD_average[:end] - D.limits[1, :end], alpha=0.7)

ax[0].set_ylabel('MSD ($nm^2$)', fontsize=14)
ax[0].set_xlabel('time (ns)', fontsize=14)
ax[0].xaxis.set_tick_params(labelsize=14)
ax[0].yaxis.set_tick_params(labelsize=14)
#ax[0].text(0.5, -1, '(a)', fontsize=14)

np.random.seed(4)  # ethanol - 4 , methanol - 5,33, 46 (16 interesting)
trajs = np.random.randint(0, D.com.shape[1], size=3)
for i in trajs:
     ax[1].plot(D.time, D.com[:, i, 2], linewidth=2)

ax[1].set_ylabel('$z$ position ($nm$)', fontsize=14)
ax[1].set_xlabel('time (ns)', fontsize=14)
ax[1].xaxis.set_tick_params(labelsize=14)
ax[1].yaxis.set_tick_params(labelsize=14)

ax[2].plot(2*D.time[:-1], D.acov, linewidth=2)

ax[2].set_ylabel('Autocovariance', fontsize=14)
ax[2].set_xlabel('k', fontsize=14)
ax[2].xaxis.set_tick_params(labelsize=14)
ax[2].yaxis.set_tick_params(labelsize=14)
ax[2].set_xlim(-1, 10)
#ax[2].set_xticklabels([-1, 1, 3, 5, 7, 9])

plt.tight_layout()
plt.savefig('msd_ztrace_acf.pdf')
plt.show()
