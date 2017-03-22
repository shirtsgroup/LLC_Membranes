#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import Thickness


def initialize():

    parser = argparse.ArgumentParser(description='Remove water layer between membrane periodic images')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file to remove water layer from')
    parser.add_argument('-o', '--output', default='wiggle_nolayer.gro', type=str, help='Name of output .gro file')
    parser.add_argument('-b', '--buffer', default=0.1, type=float, help='Float b/w 0 and 1. Used to widen the bounds on'
                                                                        'the water layer by making zmax and zmin smaller'
                                                                        'and bigger respectively')
    parser.add_argument('--index', help='If this is specified, this script will create an index group containing all of '
                                        'the waters left after they are removed', action="store_true")

    args = parser.parse_args()

    return args


def water_indices(pos, max, min, buffer):
    """
    :param pos: xyz coordinates of all atoms in system
    :param max: maximum z coordinate before water will be removed
    :param min: minimum z coordinate for water to not be removed
    :param buffer: Percent of membrane thickness to adjust max and min by
    :return: indices of water molecules to be removed
    """

    indices = [a.index for a in t.topology.atoms if a.name == 'O' and 'HOH' in str(a.residue)]

    thickness = max - min
    b = thickness * buffer
    max -= b
    min += b

    remove = []

    for i in range(len(indices)):
        if min > pos[0, indices[i], 2] or pos[0, indices[i], 2] > max:
            remove.append(indices[i])

    for i in range(len(remove)):
        remove.append(remove[i] + 1)
        remove.append(remove[i] + 2)

    return remove


def write_gro(nolayer, out):

    pos = nolayer.xyz

    ids = [a.name for a in nolayer.topology.atoms]
    res = []
    for i in range(nolayer.xyz.shape[1]):
        res.append(nolayer.topology.atom(i).residue.name)

    box = nolayer.unitcell_vectors

    atoms = pos.shape[1]
    count = 0
    with open(out, 'w') as f:
        f.write('Solvated system with water layer removed\n')
        f.write('%s\n' % atoms)
        for i in range(atoms):
            f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(count, res[i], ids[i],
                        count, pos[0, i, 0], pos[0, i, 1], pos[0, i, 2]))
            count += 1
        f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}\n'.format(box[0, 0, 0], box[0, 1, 1], box[0, 2, 2]
                                                                    ,box[0, 0, 1], box[0, 2, 0], box[0, 1, 0]
                                                                    ,box[0, 0, 2], box[0, 1, 2], box[0, 2, 0]))


def write_ndx(keep, t):
    """ Generate index groups for waters inside membrane. The indices are the same as those in the fully solvated
    structure """

    waters = []
    for a in t.topology.atoms:
        if a.index in keep and 'HOH' in str(a.residue):
            waters.append(a.index)

    print len(waters)

    count = 0
    with open('water_index.ndx', 'w') as f:
        for index in waters:
            if count % 8 != 0:
                f.write('{:<10s}'.format(str(index)))
            else:
                f.write('{:<10s}\n'.format(str(index)))
            count += 1



    return waters


if __name__ == "__main__":

    args = initialize()

    thickness, zmax, zmin = Thickness.thickness(args.gro)

    t = md.load(args.gro)

    pos = t.xyz

    remove = water_indices(pos, zmax, zmin, args.buffer)

    keep = [a.index for a in t.topology.atoms if a.index not in remove]

    nolayer = t.atom_slice(keep)

    write_gro(nolayer, args.output)

    if args.index:
        write_ndx(keep, t)
