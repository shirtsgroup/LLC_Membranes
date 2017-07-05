#! /usr/bin/env python

from builtins import range
import numpy as np
import argparse


def initialize():

    parser = argparse.ArgumentParser(description='Rename residues based on input')

    parser.add_argument('-i', '--index', default='PPPP.ndx', type=str, help='Index file describing names of residues '
                                                                            'and atomic indices included in the residue')
    parser.add_argument('-t', '--top', default='PPPP.top', type=str, help='Gromacs topology file with an '
                                                                          '[ atoms ] section')
    parser.add_argument('-ot', '--output_top', default='out.top', help='Name of output topology')
    parser.add_argument('-og', '--output_gro', default='out.gro', help='Name of output .gro coordinate file')
    parser.add_argument('-n', '--name', default='PPPP', type=str, help='Name of residue being replaced')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    with open(args.index) as f:

        residues = []
        index = []
        for line in f:
            residues.append(line.split())

    with open(args.top, 'r') as f:

        a = []
        for line in f:
            a.append(line)

    atoms_section = 0

    while a[atoms_section].count('[ atoms ]') == 0:
        atoms_section += 1

    atoms_section += 2

    res = []
    atom_no = 1
    renumber_atoms = {}
    renumber_res = {}

    for r in residues:

        if r[0] not in res:
            res.append(r[0])

        if len(r) % 2 == 1:  # a range of atomic indices are given
            for i in range(atoms_section + int(r[1]), atoms_section + int(r[2]) + 1):
                a[i] = a[i].replace(args.name, r[0])
                renumber_atoms[i - atoms_section + 1] = atom_no
                renumber_res[i - atoms_section] = len(res)
                atom_no += 1

        if len(r) % 2 == 0:  # a single atomic index is given
            index = atoms_section + int(r[1])
            a[index] = a[index].replace(args.name, r[0])
            renumber_atoms[index - atoms_section + 1] = atom_no
            renumber_res[i - atoms_section] = len(res)
            atom_no += 1

        print(len(res))

    with open(args.output_top, 'w') as f:

        for line in a:

            f.write(line)
