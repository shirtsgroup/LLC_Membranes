#!/usr/bin/env python

from LLC_Membranes.analysis.rdf import System
from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt

sol = 'URE'
pickles = ['rdf_%s.pl' % sol, 'rdf_HII_CC1C2C3C4C5.pl']
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt" % sol
opacity=0.6

print('%s/%s' % (path, pickles[0]))
rdf = file_rw.load_object('%s/%s' %(path, pickles[0]))
mean = rdf.density.mean(axis=0)
maximum = np.amax(mean[np.argwhere(rdf.r > 0.4)])
plt.plot(rdf.r, mean, lw=2)
plt.fill_between(rdf.r, rdf.errorbars[1, :] + mean, mean - rdf.errorbars[0, :], alpha=opacity)

rdf = file_rw.load_object('%s/%s' %(path, pickles[1]))
mean = rdf.density.mean(axis=0)
rdf.errorbars *= (maximum / np.max(mean))
mean *= (maximum / np.max(mean))
plt.plot(rdf.r, mean, lw=2)
plt.fill_between(rdf.r, rdf.errorbars[1, :] + mean, mean - rdf.errorbars[0, :], alpha=opacity)


plt.show()
