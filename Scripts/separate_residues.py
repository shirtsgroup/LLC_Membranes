#! /usr/bin/env python

import numpy as np
import argparse


def initialize():

    parser = argparse.ArgumentParser(description='Rename residues based on input')

    parser.add_argument('-i', '--index', default='PZPZ.ndx', type=str, help='Index file describing names of residues '
                                                                        'and atomic indices included in the residue')
    parser.add_argument('-t', '--top', help='Gromacs topology file with an [ atoms ] section')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    with open(args.index) as f:

        residues = []
        index = []
        for line in f:
            residues.append(line.split())

    with open('PZPZ.top', 'r') as f:

        a = []
        for line in f:
            a.append(line)

    atoms_section = 0

    while a[atoms_section].count('[ atoms ]') == 0:
        atoms_section += 1

    atoms_section += 2

    for r in residues:

        if len(r) % 2 == 1:  # a range of atomic indices are given
            for i in range(atoms_section + int(r[1]), atoms_section + int(r[2]) + 1):
                a[i] = a[i].replace('PZPZ', r[0])

        if len(r) % 2 == 0: # a single atomic index is given
            index = atoms_section + int(r[1])
            a[index] = a[index].replace('PZPZ', r[0])

    with open('out.top', 'w') as f:

        for line in a:

            f.write(line)
