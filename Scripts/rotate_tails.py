#! /usr/bin/env python

"""
Rotate monomer tails so they form a specified angle with the xy plane
"""

import argparse
import numpy as np
import mdtraj as md
import reposition
import tilt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the tilt angle of alkyl tails')

    parser.add_argument('-m', '--mon', default='tilted.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-b', '--bond_axes', default=["C4-O", "C3-O1", "C2-O2"], help='Bonds around which to rotate')
    parser.add_argument('-a', '--angle', default=40, type=float, help='Trajectory file (.xtc or .trr)')
    parser.add_argument('--gen_tails', help='Convenience flag for generating an index group containing the tail'
                                            'indices. This only works with the monomer which this program was written'
                                            'for', action="store_true")
    parser.add_argument('-x', '--index', type = str, default='index.ndx', help='Index file containing groups for each '
                        'tail. Each group should have a header line of the format [ groupname ]. The following line '
                        'should list the atoms in the order of their connectivity Each group should be separated by a '
                                                                               'blank line')

    args = parser.parse_args()

    return args


def atoms(bond_list, top):
    """
    :param bonds: A list of names of atoms in bonded pairs
    :return: A list with just atom names
    """

    nbonds = len(bond_list)  # number of bonds

    atoms = []  # reformat input from ["C-C"] to ["C", "C"]
    for i in range(nbonds):
        a = bond_list[i].split('-')
        for j in range(len(a)):
            atoms.append(a[j])

    indices = np.zeros([nbonds*2], dtype=int)
    for atom in top:
        if atom.name in atoms:
            ndx = atoms.index(atom.name)
            indices[ndx] = atom.index

    indices = np.reshape(indices, (nbonds, 2))

    return indices


def tail_rotate(axis, tail):

    # Step 1: Translate point on vector to origin

    t_origin = reposition.translate(tail[:, axis[0]])

    return t_origin

if __name__ == "__main__":

    args = initialize()

    full = md.load('%s' % args.mon)
    pos = full.xyz

    indices = atoms(args.bond_axes, full.topology.atoms)  # get indices of atoms in bond which we are rotating about

    tails = tilt.read_index(args.index)

    ntails = len(tails)

    for i in range(ntails):

        keep = [a.index for a in full.topology.atoms if a.name in tails[i]]
        t = full.atom_slice(keep)

        rotated = tail_rotate(indices[i, :], t.xyz)








