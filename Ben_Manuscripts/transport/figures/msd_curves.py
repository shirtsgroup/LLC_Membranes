#!/usr/bin/env python

from LLC_Membranes.analysis.msd import Diffusivity
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt

pickle = "na_msd"
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/MET/"

end = -1

fig, ax1 = plt.subplots()

inset_dimensions = [0.25, 0.6, 0.32, 0.32] # left, bottom, width, height
ax2 = fig.add_axes(inset_dimensions)

colors = ['xkcd:blue', 'xkcd:orange']

for i, wt in enumerate([10, 5]):
	
	D = file_rw.load_object('%s/%swt/%s.pl' % (path, wt, pickle))

	ax1.plot(D.time[:end], D.MSD_average[:end], color=colors[i], label='%d wt %% water' % wt)
	ax1.fill_between(D.time[:end], D.MSD_average[:end] + D.limits[0, :end], D.MSD_average[:end] - D.limits[1, :end], alpha=0.7, color=colors[i])

	if wt == 5:
		ax2.plot(D.time[:end], D.MSD_average[:end], color=colors[i])
		ax2.fill_between(D.time[:end], D.MSD_average[:end] + D.limits[0, :end], D.MSD_average[:end] - D.limits[1, :end], alpha=0.7, color=colors[i])


ax1.tick_params(labelsize=14)
ax1.set_xlabel('Time (ns)', fontsize=14)
ax1.set_ylabel('MSD (nm$^2$)', fontsize=14)
ax2.set_xlabel('Time (ns)')
ax2.set_ylabel('MSD (nm$^2$)')
leg = ax1.legend(loc=5, fontsize=13)

# adjust legend location
bb = leg.get_bbox_to_anchor().inverse_transformed(ax1.transAxes)
yoffset = -0.1
bb.y0 += yoffset
bb.y1 += yoffset
leg.set_bbox_to_anchor(bb, transform=ax1.transAxes)

plt.tight_layout()
plt.savefig('%s_comparison.pdf' % pickle)
plt.show()	
