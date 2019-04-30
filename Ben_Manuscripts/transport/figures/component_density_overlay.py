#!/usr/bin/env python

from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import numpy as np

# pickle files generated with : rdf.py -t PR_nojump.xtc -g berendsen.gro -r HII -a 'head groups' -r HII -a 'tails' -r NA -r HOH -spline

wts = [5, 10]
maxy = 0
pickles = ['NA', 'HII_tails', 'HII_head groups', 'HOH']
labels= ['Sodium', 'Tails', 'Head Groups', 'Water']
colors = ['xkcd:blue', 'xkcd:red', 'xkcd:green', 'xkcd:orange']
markers = ["s", "o"]

for j, w in enumerate(wts):
	#path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/MET/%dwt" % w
	path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/dry'

	for i, r in enumerate(pickles):

		rdf = file_rw.load_object('%s/rdf_%s.pl' % (path, r))
		mean = rdf.density.mean(axis=0)
		if j == 0:
			plt.plot(rdf.r, mean, linewidth=2, label=labels[i], color=colors[i], marker=markers[j])
		else:
			plt.plot(rdf.r, mean, linewidth=2, color=colors[i], marker=markers[j])
		plt.fill_between(rdf.r, rdf.errorbars[1, :] + mean, mean - rdf.errorbars[0, :], alpha=0.4)
		if r == 'HII_head groups':
			maximum = rdf.r[np.argmax(mean)]
			maxy = np.max(mean)
			plt.plot([maximum, maximum], [0, maxy], '--', color='black')

plt.ylabel('Density (count / nm$^3$)', fontsize=14)
plt.xlabel('Distance from pore center (nm)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=14, loc=1)
plt.tight_layout()

#plt.savefig('component_density_%swt.pdf' % wt)
plt.show()
