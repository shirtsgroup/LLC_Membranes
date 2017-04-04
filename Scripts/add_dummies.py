#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help = 'Name of coordinate file')
    parser.add_argument('-o', '--out', default='wiggled.gro', type=str, help = 'Name of output file')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)  # load trajectory using mdtraj
    pos = t.xyz  # get positions
    v = t.unitcell_vectors

    ntails = 3
    ndummies = 6
    nmonomers = 480
    Hd = ['H77', 'H78', 'H79', 'H80', 'H81', 'H82']

    atomspmon = pos.shape[1] / nmonomers

    # new_pos = np.zeros([(atomspmon + ndummies)*nmonomers, 3])

    with open(args.out, 'w') as f:

        f.write('This is a .gro file\n')
        f.write('%s\n' % ((atomspmon + ndummies)*nmonomers))

        count = 1
        count2 = 0
        for a in t.topology.atoms:
            if count2 != 0 and count2 % (atomspmon - 1) == 0 and count < (nmonomers*(atomspmon - 1 + ndummies)):
                for j in range(ndummies):
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(1, 'HII', Hd[j],
                                        count, 0, 0, 0))
                    count += 1

            f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(a.residue.index + 1, a.residue.name, a.name,
                                                    count, pos[0, count2, 0], pos[0, count2, 1], pos[0, count2, 2]))
            count += 1
            count2 += 1

        f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}\n'.format(v[0, 0, 0], v[0, 1, 1], v[0, 2, 2],
                                                                                  v[0, 0, 1], v[0, 2, 0], v[0, 1, 0],
                                                                                  v[0, 0, 2], v[0, 1, 2], v[0, 2, 0]))
