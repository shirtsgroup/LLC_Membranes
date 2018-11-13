#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Add dummy atoms so that cross-linking can occur')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-o', '--out', default='wiggled.gro', type=str, help='Name of output file')

    args = parser.parse_args()

    return args


def add_dummies(t, residue, dummy_residue, out='dummies.gro', nmon=400):

    """
    TODO: This script can be improved quite a bit, but it is functional for now while I develop xlink.py
    :param t: topology object created from mdtraj
    :return:
    """

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

    nsol = residues.count('HOH')
    natoms = t.xyz.shape[1] - nsol

    nmonomers = nmon
    Hd = ['H77', 'H78', 'H79', 'H80', 'H81', 'H82']
    ndummies = len(Hd)

    atomspmon = int(natoms / nmonomers)
    v = t.unitcell_vectors

    with open(out, 'w') as f:

        f.write('This is a .gro file\n')
        f.write('%s\n' % (int((atomspmon + ndummies) * nmonomers) + nsol))

        count = 1
        count2 = 0
        for a in t.topology.atoms:
            if count2 != 0 and count2 % (atomspmon - 1) == 0 and count < (nmonomers * (atomspmon - 1 + ndummies)):
                for j in range(ndummies):
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(int(count / atomspmon), dummy_residue,
                                                                                  Hd[j], count, 0, 0, 0))
                    count += 1

            f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(a.residue.index + 1, a.residue.name, a.name,
                                                                    count, t.xyz[0, count2, 0], t.xyz[0, count2, 1],
                                                                    t.xyz[0, count2, 2]))
            count += 1
            count2 += 1

        f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}\n'.format(v[0, 0, 0], v[0, 1, 1], v[0, 2, 2],
                                                                                  v[0, 0, 1], v[0, 2, 0], v[0, 1, 0],
                                                                                  v[0, 0, 2], v[0, 1, 2], v[0, 2, 0]))


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)  # load trajectory using mdtraj
    pos = t.xyz  # get positions
    v = t.unitcell_vectors
    add_dummies(t, 'HII', 'HIId', out=args.out)