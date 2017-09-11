#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
from llclib import transform


def initialize():

    parser = argparse.ArgumentParser(description='Rename residues based on input')

    parser.add_argument('-i', '--index', default='tail_index.ndx', type=str, help='Index file describing names of residues '
                                                                            'and atomic indices included in the residue')
    parser.add_argument('-g', '--gro', default='Halanine11V.gro', type=str, help='Name of gromacs coordinate file')
    parser.add_argument('-o', '--output', default='out.gro', help='Name of output .gro coordinate file')

    args = parser.parse_args()

    return args


class Monomer:

    def __init__(self, gro):

        t = md.load(gro)
        self.xyz = t.xyz[0]
        self.names = [a.name for a in t.topology.atoms]

        with open(gro, 'r') as f:

            a = []
            for line in f:
                a.append(line)

        P = []
        L = []
        self.res = []
        self.name = []

        # Lists of Plane atoms, line atoms, and dummy atoms
        for i in range(2, len(a) - 1):
            self.res.append(str.strip(a[i][5:10]))
            self.name.append(str.strip(a[i][10:15]))
            if a[i].count(';') != 0:
                fields = a[i].split(';')
                annotations = fields[1].split()
                if 'P' in annotations:
                    P.append(int(i - 2))
                if 'L' in annotations:
                    L.append(int(i - 2))  # adjust for top lines and index (count from 0 rather than 1)
                if 'R' in annotations:
                    self.ref_atom_index = int(i - 2)

        self.plane = np.zeros([3, 3])

        for i in range(len(P)):
            self.plane[i, :] = self.xyz[P[i], :]

    def writeout(self, outname):
        """
        :param outname: Name of output file
        :return: write an output .gro file from coordinates
        """

        with open(outname, 'w') as f:

            f.write('This is a .gro file\n')
            f.write('%s\n' % self.xyz.shape[0])

            for i in range(self.xyz.shape[0]):

                f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format((i + 1) % 100000, '%s' % self.res[i],
                                '%s' % self.name[i], (i + 1) % 100000, self.xyz[i, 0], self.xyz[i, 1], self.xyz[i, 2]))

            f.write('{:10f}{:10f}{:10f}\n'.format(0, 0, 0))

    def tails(self, indices):


if __name__ == "__main__":

    args = initialize()

    mon = Monomer(args.gro)
    mon.xyz = transform.rotateplane_coords(mon.xyz, mon.plane, angle=0)
    mon.tails(args.index)
    mon.writeout(args.output)