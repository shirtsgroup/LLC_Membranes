#!/usr/bin/env python

import mdtraj as md
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import tqdm

wt = 10 
water = True 
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/pure_water/%swt" % wt

t = md.load('%s/PR_box.xtc' % path, top='%s/berendsen.gro' % path)
print('Trajectory Loaded')

xmin = t.xyz[..., 0].min()
xmax = t.xyz[..., 0].max()
ymin = t.xyz[..., 1].min()
ymax = t.xyz[..., 1].max()

rang = [[xmin, xmax], [ymin, ymax]]

if water:
	w = [a.index for a in t.topology.atoms if a.residue.name == 'HOH']
	pos = t.xyz[:, w, :]
else:
	pos = t.xyz

bins = 100
hist = np.zeros([bins, bins])
for f in tqdm.tqdm(range(t.n_frames)):
	H, xedges, yedges = np.histogram2d(pos[f, :, 0], pos[f, :, 1], bins=bins, range=rang)
	hist += H

hist /= t.n_frames

ax = plt.gca()
cmap = 'viridis'
im = ax.imshow(hist.T, extent=[xmin, xmax, ymin, ymax], interpolation='gaussian', cmap=cmap, vmin=0, vmax=hist[np.nonzero(hist)].max())
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im, cax=cax)

if water:
	cbar.set_label('Density ($H_2O / nm^2$)', fontsize=14)
else:
	cbar.set_label('Density ($Atoms / nm^2$)', fontsize=14)

ax.set_xlabel('$x$-location (nm)', fontsize=14)
ax.set_ylabel('$y$-location (nm)', fontsize=14)
ax.tick_params(labelsize=14)
plt.tight_layout()
if water:
	plt.savefig('total_water_density_%swt.pdf' % wt)
else:
	plt.savefig('total_density_%swt.pdf' % wt)
plt.show()
