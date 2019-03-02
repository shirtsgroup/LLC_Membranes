#!/usr/bin/env python

from LLC_Membranes.analysis.ztrace import ZTrace
from LLC_Membranes.llclib import file_rw
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
import numpy as np

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/THF/10wt/"
THF_resnumber = 14
NA_resnumber = 327

THF = file_rw.load_object('%s/trace.pl' % path)
NA = file_rw.load_object('%s/na_trace.pl' % path)

THF_trace = THF.com[:, THF_resnumber, 2]
THF_rd = THF.radial_distance[:, THF_resnumber]

NA_trace = NA.com[:, NA_resnumber, 2]
NA_rd = NA.radial_distance[:, NA_resnumber]

fig, ax = plt.subplots(2, 1, sharex=True, figsize=(9, 9))

colormap = 'plasma_r'
norm = plt.Normalize(0, 1.5)

# THF trace
points = np.array([THF.t.time / 1000, THF_trace]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = LineCollection(segments, cmap=colormap, norm=norm)
lc.set_array(THF_rd)
lc.set_linewidth(2)
line1 = ax[0].add_collection(lc)

# Sodium trace
points = np.array([NA.t.time / 1000, NA_trace]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lc = LineCollection(segments, cmap=colormap, norm=norm)
lc.set_array(NA_rd)
lc.set_linewidth(2)
line2 = ax[0].add_collection(lc)

ax[0].set_xlim(0, THF.t.time[-1] / 1000)
ymin = min(np.amin(THF_trace), np.amin(NA_trace))
ymax = max(np.amax(THF_trace), np.amax(NA_trace))
span = ymax - ymin
ymax += .05 * span
ymin -= .05 * span
ax[0].set_ylim(ymin, ymax)

ax[0].set_ylabel('COM $z$-position (nm)', fontsize=14)
ax[0].tick_params(labelsize=14)

divider = make_axes_locatable(ax[0])
cax = divider.append_axes("top", size="7%", pad="15%")
#cax.xaxis.set_ticks_position("top")
cb = colorbar(line1, cax=cax, orientation="horizontal")
#cbar = fig.colorbar(line, ax=cax)#, orientation="horizontal", location='top')
#cb.set_label('Radial distance from pore center (nm)', fontsize=14)

t = md.load('%s/PR.trr' % path, top='%s/berendsen.gro' % path)

thf_index = [a.index for a in t.topology.atoms if a.residue.name == 'THF' and a.name == 'O'][THF_resnumber]

NA_index = [a.index for a in t.topology.atoms if a.residue.name == 'NA'][NA_resnumber]


separation = np.linalg.norm(t.xyz[:, thf_index, :] - t.xyz[:, NA_index, :], axis=1)

ax[1].plot(t.time / 1000, separation, linewidth=2)
ax[1].set_ylim(0, np.amax(separation))
ax[1].set_xlabel('Time (ns)', fontsize=14)
ax[1].tick_params(labelsize=14)
ax[1].set_ylabel('Sodium-THF distance', fontsize=14)

plt.tight_layout()
plt.savefig('thf_sodium_coordination.pdf') 
plt.show()

