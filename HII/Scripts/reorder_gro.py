#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Rewrite .gro file so that it matches topology')  # allow input from user

    parser.add_argument('-i', '--input', default='HII_packed.gro', type=str, help='Name of input file')
    parser.add_argument('-s', '--solvent', default='SOL', help='Name of solvent')
    parser.add_argument('-o', '--output', default='Ordered.gro', help='Name of reordered output file')
    parser.add_argument('-a', '--natoms', default=138, type=int, help='Number of atoms per molecule')
    parser.add_argument('-I', '--ion', default='NA', type=str, help='Name of ion in system')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.input)

    # The following works, looks cleaner, but is annoyingly slow
    # ions = [a.index for a in t.topology.atoms if a.name != args.ion]
    # residue = [a for a in range(t.n_atoms) if a not in ions]

    with open(args.input, 'r') as f:

        a = []
        ions = 0
        for line in f:
            a.append(line)
            if line.count(args.ion) > 0:
                ions += 1

    nres = int((len(a) - 3) / args.natoms)
    ions_per = int(ions / nres)
    nres_atoms = args.natoms - ions_per

    f = open(args.output, 'w')

    f.write('This is a .gro file\n')
    f.write('%s\n' % (len(a) - 3))

    line = 2
    count = 1
    for i in range(nres):
        for j in range(nres_atoms):
            f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(i + 1, str(a[line][5:10]),
            str(a[line][10:15]), count, float(a[line][20:28]), float(a[line][28:36]), float(a[line][36:44])))
            line += 1
            count += 1
        line += ions_per

    line = 2 + nres_atoms
    for i in range(nres):
        for j in range(ions_per):
            f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(nres + i*ions_per + j + 1, str(a[line][5:10]),
            str(a[line][10:15]), count, float(a[line][20:28]), float(a[line][28:36]), float(a[line][36:44])))
            line += 1
            count += 1
        line += nres_atoms

    f.write(a[-1])  # box dimensions

    f.close()
