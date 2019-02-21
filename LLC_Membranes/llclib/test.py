#!/usr/bin/env python

from matplotlib import path
import numpy as np
import mdtraj as md

# parameters from gaff.dat

sigma = 1.4870
epsilon = 0.0157

# convert to gromacs units

print 'Sigma: %1.6f' % (2*sigma/(10*2**(1.0/6)))
print 'Epsilon: %1.6f' %(4.184 * epsilon)

