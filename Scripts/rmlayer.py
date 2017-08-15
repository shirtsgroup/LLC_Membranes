#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import Thickness
from llclib import top, file_rw, physical
from hydrophobic_density import limits


def initialize():

    parser = argparse.ArgumentParser(description='Remove water layer between membrane periodic images')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file to remove water layer from')
    parser.add_argument('-o', '--output', default='wiggle_nolayer.gro', type=str, help='Name of output .gro file')
    parser.add_argument('-b', '--buffer', default=0.1, type=float, help='Float b/w 0 and 1. Used to widen the bounds on'
                                                                        'the water layer by making zmax and zmin smaller'
                                                                        'and bigger respectively')
    parser.add_argument('--index', help='If this is specified, this script will create an index group containing all of '
                                        'the waters left after they are removed', action="store_true")
    parser.add_argument('-t', '--tail', default=['O5', 'O6', 'O7', 'O8', 'O9', 'O10'], help='Name of atoms in tail to'
                                                'be used as a reference for waters that should be kept in the region')
    parser.add_argument('-p', '--pore', default=['C1', 'C2', 'C3', 'C4', 'C5', 'C'], help='Name of atoms in pore region'
                                                'to be used as reference for waters that should stay in pore. Dont use'
                                                'an ion for this since there are ions in solution')

    args = parser.parse_args()

    return args


def water_indices(pos, max, min, buffer, t):
    """
    :param pos: xyz coordinates of all atoms in system
    :param max: maximum z coordinate before water will be removed
    :param min: minimum z coordinate for water to not be removed
    :param buffer: Percent of membrane thickness to adjust max and min by
    :param t: mdtraj topology object
    :return: indices of water molecules to be removed
    """

    indices = [a.index for a in t.topology.atoms if a.name == 'O' and 'HOH' in str(a.residue)]
    #indices = [a.index for a in t.topology.atoms if a.residue.is_water]
    indices = t.topology.select("water")

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
        # map the names back to the tip3p names...need to figure out a better way
        if res[i] == 'HOH':
            res[i] = 'SOL'
            if ids[i] == 'H1':
                ids[i] = 'HW1'
            if ids[i] == 'H2':
                ids[i] = 'HW2'
            if ids[i] == 'O':
                ids[i] = 'OW'

    box = nolayer.unitcell_vectors

    atoms = pos.shape[1]
    count = 1
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
    membrane = []
    for a in t.topology.atoms:
        if a.index in keep and 'HOH' in str(a.residue):  # if the atom is being kept and is part of water, record it
            waters.append(a.index)
        elif a.index in keep:  # otherwise it is part of the membrane. Needs to be in keep though or else the unkept \
            membrane.append(a.index)  # water will go in the membrane list where they aren't supposed to >:(

    count = 1
    with open('water_index.ndx', 'w') as f:  # open up an index file to write to

        f.write('[  water  ]\n')  # first index group
        for index in waters:
            if count % 10 != 0:  # every 10 entries, make a new line
                f.write('{:<8s}'.format(str(index + 1)))  # things are indexed starting at 0 in mdtraj and 1 in gromacs
            else:
                f.write('{:<8s}\n'.format(str(index + 1)))
            count += 1

        f.write('\n[  membrane  ]\n')  # membrane section!
        count = 1
        for index in membrane:
            if count % 10 != 0:
                f.write('{:<8s}'.format(str(index + 1)))
            else:
                f.write('{:<8s}\n'.format(str(index + 1)))
            count += 1


if __name__ == "__main__":

    args = initialize()

    thickness, zmax, zmin = physical.thickness(args.gro, ref_atoms=args.tail)

    t = md.load(args.gro)

    # create separate objects for the different groups of atoms we are interested in
    pos = t.xyz  # all positions
    # get indices of all atoms of interest
    water_index = t.topology.select("water")
    pore_index = [a.index for a in t.topology.atoms if a.name in args.pore]
    tail_index = [a.index for a in t.topology.atoms if a.name in args.tail]

    # create new objects by slicing from the full trajectory object and only keep their coordinates (.xyz)
    water = t.atom_slice(water_index).xyz
    pore = t.atom_slice(pore_index).xyz
    tail = t.atom_slice(tail_index).xyz

    p_centers = physical.avg_pore_loc(4, pore)  # calculate pore centers based on average sodium ion positions
    # print p_centers

    pore_limits, pore_std = limits(pore, p_centers)
    tail_limits, tail_std = limits(tail, p_centers)
    # print pore_limits
    # print tail_limits

    remove = water_indices(pos, zmax, zmin, args.buffer, t)

    # keep = [a.index for a in t.topology.atoms if a.index not in remove]
    #
    # nolayer = t.atom_slice(keep)
    nolayer, keep = top.keep(remove, t, exclude=True)

    write_gro(nolayer, args.output)

    if args.index:
        file_rw.write_water_ndx(keep, t)

