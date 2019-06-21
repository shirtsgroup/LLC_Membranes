#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
from LLC_Membranes.llclib import topology


def initialize():

    parser = argparse.ArgumentParser(description='Add dummy atoms so that cross-linking can occur')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-o', '--out', default='wiggled.gro', type=str, help='Name of output file')
    parser.add_argument('-r', '--residue', default='HII', type=str, help='Name of residue to be replaced by a dummy'
                                                                         'residue')
    parser.add_argument('-d', '--dummy_residue', default='HIId', type=str, help='Name of dummy residue')

    args = parser.parse_args()

    return args


def add_dummies(gro, residue, dummy_residue, out='dummies.gro'):
    # TODO: This script can be improved quite a bit, but it is functional for now while I develop xlink.py
    """ Add dummy hydrogen atoms to chosen residues in a configuration

    :param t: topology object created from mdtraj
    :param residue: name of residue to which dummy atoms are added
    :param dummy_residue: name of dummy_residue associated with dummy topology file
    :param out: name of output .gro file

    :type t: mdtraj object
    :type residue: str
    :type dummy_residue: str
    :type out: str

    :return: .gro file with dummies atoms added. Note: an energy minimization is necessary to snap the dummies into
    place
    """

    t = md.load(gro)
    names = topology.fix_names(gro)
    for i, a in enumerate(t.topology.atoms):
        a.name = names[i]

    original_residue = topology.LC('%s' % residue)
    LC = topology.LC('%s' % dummy_residue)

    residues = [a.residue.name for a in t.topology.atoms]

    if 'HOH' in residues:
        for a in t.topology.atoms:
            if a.residue.name == 'HOH':
                if a.name == 'O':
                    a.name = 'OW'
                elif a.name == 'H1':
                    a.name = 'HW1'
                elif a.name == 'H2':
                    a.name = 'HW2'
        for a in t.topology.atoms:  # if you change residue name before changing names, this won't work
            if a.residue.name == 'HOH':
                a.residue.name = 'SOL'
            if a.residue.name == residue:
                a.residue.name = dummy_residue

    nsol = 0
    for i in set(residues):
        if i != residue:
            nsol += residues.count(i)

    natoms = t.n_atoms - nsol

    Hd = LC.dummies
    ndummies = len(Hd)

    atomspmon = original_residue.natoms
    nmonomers = natoms // atomspmon
    v = t.unitcell_vectors

    with open(out, 'w') as f:

        f.write('This is a .gro file\n')
        f.write('%s\n' % (int((atomspmon + ndummies) * nmonomers) + nsol))

        count = 1
        count2 = 0
        for a in t.topology.atoms:
            if count2 != 0 and count2 % atomspmon == 0 and count < (nmonomers * (atomspmon + ndummies)):
                for j in range(ndummies):
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(int(count2 / atomspmon), dummy_residue,
                                                                                  Hd[j], count, 0, 0, 0))
                    count += 1

            # kind of hacky
            if a.residue.name == original_residue.residues[0]:
                a.residue.name = dummy_residue

            f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(a.residue.index + 1, a.residue.name, a.name,
                                                                          count % 100000, t.xyz[0, count2, 0],
                                                                          t.xyz[0, count2, 1], t.xyz[0, count2, 2]))
            count += 1
            count2 += 1

        f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}\n'.format(v[0, 0, 0], v[0, 1, 1], v[0, 2, 2],
                                                                                  v[0, 0, 1], v[0, 2, 0], v[0, 1, 0],
                                                                                  v[0, 0, 2], v[0, 1, 2], v[0, 2, 0]))


if __name__ == "__main__":

    args = initialize()

    add_dummies(args.gro, args.residue, args.dummy_residue, out=args.out)
