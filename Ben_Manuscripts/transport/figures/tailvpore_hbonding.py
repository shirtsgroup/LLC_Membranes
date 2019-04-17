#!/usr/bin/env python

import numpy as np
from LLC_Membranes.analysis import hbonds
from LLC_Membranes.llclib import topology, physical

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/pure_water"
spline = True
tails = True  # look at hbonds in the tails
wt = 5
r = 1.5  # cut-off between distal tails
npores = 4
dwell_fraction = 0.95

full_path = path + '/%swt' % wt

hb = hbonds.System('%s/PR_nojump.xtc' % full_path, '%s/berendsen.gro' % full_path)

hb.set_eligible('HOH', 'all')

# find pore centers
pore_defining_atoms = topology.LC('NAcarb11V').pore_defining_atoms
pore_atoms = [a.index for a in hb.t.topology.atoms if a.name in pore_defining_atoms]

if spline:
    print('Creating pore splines')
    pore_centers = physical.trace_pores(hb.t.xyz[:, pore_atoms, :], hb.t.unitcell_vectors, 10)[0]
else:
    pore_centers = physical.avg_pore_loc(npores, hb.t.xyz[:, pore_atoms, :], hb.t.unitcell_vectors)[0]

oxygen = [a.index for a in hb.t.topology.atoms if a.residue.name == 'HOH' and a.name == 'OW']
inregion = physical.partition(hb.t.xyz[:, oxygen, :], pore_centers, r, buffer=0, unitcell=hb.t.unitcell_vectors, npores=npores, spline=True)

if tails:
    inregion = ~inregion  # '~' flips True and False

dwell = np.full((hb.t.n_frames, len(oxygen)), False, dtype=bool)

for t in range(hb.t.n_frames):
    dwell[t, inregion[t]] = True

fraction_dwelled = np.sum(dwell, axis=0) / hb.t.n_frames  # fraction of total time spend in region of interest

keep = np.where(fraction_dwelled >= dwell_fraction)[0]

# oxygen, hb.A, and hb.D are ordered the same way
hb.A = np.array(hb.A)[keep]
hb.D = np.array(hb.D)[keep]
from itertools import chain
hkeep = list(chain.from_iterable((2*i, 2*i + 1) for i in keep))  # since there are two H's per O
hb.H = np.array(hb.H)[hkeep]

total_waters = len(hb.A)

hb.identify_hbonds(0.35, 30)
n = [a.shape[1] for a in hb.hbonds]

print('Fraction of water molecules involved in an hbond: %s' % (np.mean(n) / total_waters))
