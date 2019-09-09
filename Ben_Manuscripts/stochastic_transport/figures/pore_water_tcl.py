#!/usr/bin/env python

import mdtraj as md
import numpy as np

from LLC_Membranes.llclib import physical, topology

r = 1

t = md.load('initial.gro')
keep = [a.index for a in t.topology.atoms if a.residue.name == 'HOH']
res_start = keep[0]


com = physical.center_of_mass(t.xyz[:, keep, :], [18., 1., 1.])
membrane = topology.LC('HII')  # object w/ attributes of LC making up membrane
hg = [a.index for a in t.topology.atoms if a.name in membrane.pore_defining_atoms and a.residue.name
        == membrane.name]
pore_centers = physical.avg_pore_loc(4, t.xyz[:, hg, :], t.unitcell_vectors, buffer=0, spline=False)

partition = physical.partition(com, pore_centers, r, unitcell=t.unitcell_vectors,
                                            spline=False)

pore_indices = [res_start + 3 * i for i in np.where(partition[0, :])[0]]
tail_indices = [res_start + 3 * i for i in np.where(partition[0, :] == False)[0]]  # have to use double equals sign. Using is doesn't work with np.where

with open('partition.tcl', 'w') as f:
    f.write('color Display Background white\n')
    f.write('mol addrep 0\n')
    f.write('mol modselect 0 0 index')
    for i in pore_indices:
        end = i + 3
        f.write(' %s to %s' % (i, end - 1))
    f.write('\n')
    f.write('mol modcolor 0 0 ColorID 0\n')
    f.write('mol modstyle 0 0 CPK 2.0 0.3 12.0 12.0\n')
    f.write('mol addrep 0\n')
    f.write('mol modselect 1 0 index')
    for i in tail_indices:
        end = i + 3
        f.write(' %s to %s' % (i, end - 1))
    f.write('\n')
    f.write('mol modstyle 1 0 CPK 2.0 0.3 12.0 12.0\n')
    f.write('mol modcolor 1 0 ColorID 1\n')

