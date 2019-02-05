#!/usr/bin/env python

import numpy as np
from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import names

r = 'BUT'
atoms = ['O', 'C']
wt=10
maximum = 0


path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" %(r,wt)

for a in atoms:

	rdf = file_rw.load_object('%s/rdf_%s_%s.pl' %(path, r, a))

	mean = rdf.density.mean(axis=0)
	max_loc = np.argmax(mean)
	plt.plot([rdf.r[max_loc], rdf.r[max_loc]], [0, 1.05*np.max(mean)], '--', color='black', linewidth=2)
	plt.plot(rdf.r, mean, label='%s : %s' % (names.res_to_name[r], a))
	plt.fill_between(rdf.r, rdf.errorbars[1, :] + mean, mean - rdf.errorbars[0, :], alpha=0.7)

plt.ylabel('Density (count / nm$^3$)', fontsize=14)
plt.xlabel('Distance from pore center (nm)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('butanol_CO.pdf')
plt.show()

