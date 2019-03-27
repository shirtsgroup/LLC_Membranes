#!/usr/bin/env python

from LLC_Membranes.analysis.msd import Diffusivity
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt

# pickles generated with following command in pure_water directories
# msd.py -t PR_nojump.xtc -g berendsen.gro -r NA/water -s na/water_msd -F 0.4 --nofit -pores

pickle = "na_msd"
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/pure_water/"

end = 800 

fig, ax1 = plt.subplots()

inset_dimensions = [0.225, 0.6, 0.32, 0.32] # left, bottom, width, height
inset_dimensions = [0.24, 0.64, 0.3, 0.3]  # for sodium
ax2 = fig.add_axes(inset_dimensions)

colors = ['xkcd:blue', 'xkcd:orange']
opacity = 0.4

for i, wt in enumerate([10, 5]):
	
	D = file_rw.load_object('%s/%swt/%s.pl' % (path, wt, pickle))

	ax1.plot(D.time[:end], D.MSD_average[:end], color=colors[i], label='%d wt %% water' % wt)
	ax1.fill_between(D.time[:end], D.MSD_average[:end] + D.limits[0, :end], D.MSD_average[:end] - D.limits[1, :end], alpha=opacity, color=colors[i])

	if wt == 5:
		ax2.plot(D.time[:end], D.MSD_average[:end], color=colors[i])
		ax2.fill_between(D.time[:end], D.MSD_average[:end] + D.limits[0, :end], D.MSD_average[:end] - D.limits[1, :end], alpha=opacity, color=colors[i])


ax1.tick_params(labelsize=14)
ax1.set_xlabel('Time (ns)', fontsize=14)
ax1.set_ylabel('MSD (nm$^2$)', fontsize=14)
ax2.set_xlabel('Time (ns)')
ax2.set_ylabel('MSD (nm$^2$)')
leg = ax1.legend(loc=5, fontsize=13)

# adjust legend location
bb = leg.get_bbox_to_anchor().inverse_transformed(ax1.transAxes)
yoffset = -0.1
yoffset = -0.05  # for water
bb.y0 += yoffset
bb.y1 += yoffset
leg.set_bbox_to_anchor(bb, transform=ax1.transAxes)

plt.tight_layout()
plt.savefig('%s_comparison.pdf' % pickle)
plt.show()	
