#! /usr/bin/env python

import argparse
import mdtraj as md
from llclib import file_rw
import subprocess


def initialize():

    parser = argparse.ArgumentParser(description='Expand box and translate atoms to the z center for water gap creation')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of gro file whose box will be modified')
    parser.add_argument('-l', '--gap', default=3, type=float, help='Angle to rotate by')
    parser.add_argument('-o', '--output', default='gap.gro', type=str, help='Name of new .gro produced')
    parser.add_argument('-s', '--solvate', action="store_true", help='Solvate system after adding gap')
    parser.add_argument('-t', '--top', default='topol.top', help='name of topology to modify if solvating')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)

    box = t.unitcell_vectors
    box[0, 2, 2] += args.gap
    t.unitcell_vectors = box

    pos = t.xyz[0]
    for i in range(pos.shape[0]):
        pos[i, 2] += (args.gap / 2)

    t.xyz = pos

    file_rw.write_gro(t, args.output)

    if args.solvate:

        subprocess.call(['gmx', 'solvate', '-cp', '%s' % args.output, '-cs', 'spc216.gro', '-p', '%s' % args.top, '-o', 'solvated.gro'])

        # t = md.load('solvated.gro')
        #
        # pos = t.xyz
        #
        # water = [a.index for a in t.topology.atoms if a.residue == 'SOL']
        #
        # print(len(water))
