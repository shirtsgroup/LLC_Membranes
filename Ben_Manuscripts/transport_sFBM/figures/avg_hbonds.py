#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw, topology
from LLC_Membranes.analysis.hbonds import System, Topology, Residue
import names


def calculate_moving_average(data, n):
    """ Calculate moving average of a time series

    :param n: Number of previous points to average

    :type n: int
    """

    ret = np.cumsum(data, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def smooth(y, box_pts):
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
#solutes = ["GLY", "GCL", "PG"]
solutes = ["TET", "RIB"]
solutes = ["GCL", "SOH", "GLY", "DMP"]
nsolutes = 24
npts_ma = 8  # number of points to average over in moving average

for i in solutes:
    sys = file_rw.load_object('%s/%s/10wt/hbonds.pl' % (path, i))
    n = [a.shape[1] / nsolutes for a in sys.hbonds]
    #plt.plot(sys.t.time[:-(npts_ma - 1)] / 1000, calculate_moving_average(n, npts_ma))
    #plt.plot(sys.t.time / 1000, n)
    plt.plot(sys.t.time / 1000, smooth(n, npts_ma), label=names.res_to_name[i])

#plt.ylim(0, 1)
plt.xlabel('Time (ns)', fontsize=14)
plt.ylabel('Number of hydrogen bonds per solute', fontsize=14)
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend()
plt.tight_layout()
plt.show()
