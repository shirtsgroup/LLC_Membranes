#!/usr/bin/env python

import numpy as np
from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import names
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition

r = 'BUT'
atoms = ['C', 'O']
wt=10
maximum = 0

fig, ax = plt.subplots()

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%dwt" %(r,wt)

for a in atoms:

	rdf = file_rw.load_object('%s/rdf_%s_%s.pl' %(path, r, a))

	mean = rdf.density.mean(axis=0)
	max_loc = np.argmax(mean)
	ax.plot([rdf.r[max_loc], rdf.r[max_loc]], [0, 1.05*np.max(mean)], '--', color='black', linewidth=2)
	ax.plot(rdf.r, mean, label='%s : %s' % (names.res_to_name[r], a))
	ax.fill_between(rdf.r, rdf.errorbars[1, :] + mean, mean - rdf.errorbars[0, :], alpha=0.7)

#im = plt.imread('butanol_labeled.png')
#axins = plt.axes([0, 0, 1, 1])
#ip = InsetPosition(ax, [0.58, 0.45, 0.4, 0.4])
#axins = inset_axes(ax, width="100%", height="100%", loc='center right')
#axins.set_axes_locator(ip)
#axins.imshow(im)
#axins.axis('off')

ax.set_ylabel('Density (count / nm$^3$)', fontsize=14)
ax.set_xlabel('Distance from pore center (nm)', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
ax.legend(fontsize=14)
plt.tight_layout()
#plt.savefig('%s_%s.pdf' % (names.res_to_name[r], ''.join(atoms)))
#plt.savefig('butanol_CO.pdf')
plt.show()

