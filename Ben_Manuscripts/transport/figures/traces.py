#!/usr/bin/env python

from LLC_Membranes.analysis.ztrace import ZTrace
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt

residue = "SOH"
indices = [17, 3]  # residue number whose z trace to plot
wt = 10

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/%swt/trace.pl" % (residue, wt)

trace = file_rw.load_object(path)

trace.plot_trace(indices, cmax=1.5, savename="%s_trace.pdf" % residue)
