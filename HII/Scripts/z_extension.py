#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from llclib import file_rw


def initialize():

    parser = argparse.ArgumentParser(description='Calculate a correlation function of positions of atoms in the z direction')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-n', '--nimages', default=3, type=int, help='Number of times to duplicate unit cell in z dir')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()
    t = md.load(args.gro)
    full_box = t.unitcell_vectors
    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], args.nimages*full_box[0, 2, 2], full_box[0, 0, 1],
                   full_box[0, 2, 0], full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]
    res = [a.residue.name for a in t.topology.atoms]*args.nimages
    ids = [a.name for a in t.topology.atoms]*args.nimages
    pos = t.xyz[0, :, :]
    natoms = pos.shape[0]
    extended = np.zeros([natoms*args.nimages, 3])
    z = full_box[0, 2, 2]

    for i in range(args.nimages):
        extended[i*natoms:(i+1)*natoms, :] = pos + [0, 0, i*z]

    file_rw.write_gro_pos(extended, 'extended.gro', box=box_gromacs, res=res, ids=ids)