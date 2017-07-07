#! /usr/bin/env python

import argparse
import numpy as np
import subprocess


def initialize():

    parser = argparse.ArgumentParser(description='Create random configurations with certain features while the rest of '
                                                 'space is filled uniformly with particles')

    parser.add_argument('-i', '--input', default='wiggle.gro', type=str, help='Name of input .gro file')
    parser.add_argument('-o', '--output', default='wiggle.pdb', type=str, help='Name of output .pdb file')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    # Make gromacs do most of the hard stuff
    subprocess.call(['gmx', 'editconf', '-f', '%s' % args.input, '-o', '%s' % args.output])

    a = []
    with open(args.output, 'r') as f:
        for line in f:
            a.append(line)

    out = []
    atom_start = 0
    while a[atom_start].count('ATOM') == 0:
        out.append(a[atom_start])
        atom_start += 1

    while a[atom_start].count('TER') == 0:
        replacement = a[atom_start].split()
        atom = replacement[2]
        if atom == 'NA':
            element = 'NA'
        else:
            element = atom[0]
        replacement.append(element)
        out.append(a[atom_start].rstrip() + '{:>12}\n'.format(element))

        atom_start += 1

    out.append('TER\n')
    out.append('ENDMDL')

    with open(args.output, 'w') as f:

        for line in out:
            f.write(line)