#!/usr/bin/env python

from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import numpy as np

# pickle files generated with : rdf.py -t PR_nojump.xtc -g berendsen.gro -r HII -a 'head groups' -r HII -a 'tails' -r NA -r HOH -spline

wt = 10 
dry = True
opacity = 0.2

if dry:
	path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/dry"
	pickles = ['NA', 'HII_tails', 'HII_head groups']
else:
	path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/pure_water/%dwt" % wt
	pickles = ['NA', 'HII_tails', 'HII_head groups', 'HOH']
labels= ['Sodium', 'Tails', 'Head Groups', 'Water']
colors = ['xkcd:blue', 'xkcd:red', 'xkcd:green', 'xkcd:orange']

for i, r in enumerate(pickles):

	rdf = file_rw.load_object('%s/rdf_%s.pl' % (path, r))
	print(r)
	mean = rdf.density.mean(axis=0)
	plt.plot(rdf.r, mean, linewidth=2, label=labels[i], color=colors[i])
	plt.fill_between(rdf.r, rdf.errorbars[1, :] + mean, mean - rdf.errorbars[0, :], alpha=opacity)
	if r == 'HII_head groups':
		maximum = rdf.r[np.argmax(mean)]
		maxy = np.max(mean)

plt.ylabel('Density (count / nm$^3$)', fontsize=14)
plt.xlabel('Distance from pore center (nm)', fontsize=14)
plt.plot([maximum, maximum], [0, maxy], '--', color='black')
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=14, loc=1)
plt.tight_layout()

if dry:
	plt.savefig('component_density_dry.pdf')
else:
	plt.savefig('component_density_%swt.pdf' % wt)
plt.show()
