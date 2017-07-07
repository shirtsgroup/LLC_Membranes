#! /usr/bin/env python

# from builtins import range
from __future__ import print_function
import numpy as np
import argparse


def initialize():

    parser = argparse.ArgumentParser(description='Rename residues based on input')

    parser.add_argument('-i', '--index', default='PZPZ.ndx', type=str, help='Index file describing names of residues '
                                                                            'and atomic indices included in the residue')
    parser.add_argument('-t', '--top', default='PZPZ.top', type=str, help='Gromacs topology file with an '
                                                                          '[ atoms ] section')
    parser.add_argument('-g', '--gro', default='PZPZ.gro', type=str, help='Name of gromacs coordinate file')
    parser.add_argument('-ot', '--output_top', default='out.top', help='Name of output topology')
    parser.add_argument('-og', '--output_gro', default='out.gro', help='Name of output .gro coordinate file')
    parser.add_argument('-n', '--name', default='PZPZ', type=str, help='Name of residue being replaced')

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

    # now rewrite topology
    out = []
    atoms_section = 0

    while a[atoms_section].count('[ atoms ]') == 0:
        out.append(a[atoms_section])
        atoms_section += 1

    out.append(a[atoms_section])
    out.append(a[atoms_section + 1])
    atoms_section += 2

    res = []
    atom_no = 1
    renumber_atoms = {}
    renumber_res = {}

    for r in residues:

        if r[0] not in res:
            res.append(r[0])
            res_no = len(res)
        else:
            res_no = res.index(r[0]) + 1

        if len(r) % 2 == 1:  # a range of atomic indices are given
            for i in range(atoms_section + int(r[1]), atoms_section + int(r[2]) + 1):
                a[i] = a[i].replace(args.name, r[0])
                renumber_atoms[i - atoms_section + 1] = atom_no
                renumber_res[i - atoms_section + 1] = res_no
                atom_no += 1

        if len(r) % 2 == 0:  # a single atomic index is given
            index = atoms_section + int(r[1])
            a[index] = a[index].replace(args.name, r[0])
            renumber_atoms[index - atoms_section + 1] = atom_no
            renumber_res[index - atoms_section + 1] = res_no
            atom_no += 1

    used_residues = []
    for i in range(atoms_section, atoms_section + atom_no - 1):
        k = list(renumber_atoms.keys())[list(renumber_atoms.values()).index(i - atoms_section + 1)]
        replacement = a[k + atoms_section - 1].split()
        replacement[0] = str(i - atoms_section + 1)
        replacement[2] = str(renumber_res[k])  # residue number associated with that atom number
        replacement[5] = str(i - atoms_section + 1)

        if replacement[3] not in used_residues:
            used_residues.append(replacement[3])
            out.append('; residue %s : %s\n' % (renumber_res[k], replacement[3]))

        replacement[3] = replacement[3][0]
        out.append('{:>6}{:>5}{:>6}{:>6}{:>6}{:>6}{:>13}{:>13}\n'.format(*tuple(replacement)))

    ################# BONDS #######################

    out.append('\n[ bonds ]\n')
    bonds_section = 0

    while a[bonds_section].count('[ bonds ]') == 0:
        bonds_section += 1

    bonds_section += 1
    out.append(a[bonds_section])
    bonds_section += 1

    while a[bonds_section] != '\n':
        replacement = a[bonds_section].split()
        replacement[0] = renumber_atoms[int(replacement[0])]
        replacement[1] = renumber_atoms[int(replacement[1])]
        out.append('{:>6}{:>6}{:>6}{:>13}{:>13}\n'.format(*tuple(replacement[:5])))
        bonds_section += 1

    ################## PAIRS #####################

    out.append('\n[ pairs ]\n')
    pairs_section = 0

    while a[pairs_section].count('[ pairs ]') == 0:
        pairs_section += 1

    pairs_section += 1
    out.append(a[pairs_section])
    pairs_section += 1

    while a[pairs_section] != '\n':
        replacement = a[pairs_section].split()
        replacement[0] = renumber_atoms[int(replacement[0])]
        replacement[1] = renumber_atoms[int(replacement[1])]
        out.append('{:>6}{:>6}{:>6}\n'.format(*tuple(replacement[:3])))
        pairs_section += 1

    ################## ANGLES #####################

    out.append('\n[ angles ]\n')
    angles_section = 0

    while a[angles_section].count('[ angles ]') == 0:
        angles_section += 1

    angles_section += 1
    out.append(a[angles_section])
    angles_section += 1

    while a[angles_section] != '\n':
        replacement = a[angles_section].split()
        for i in range(3):
            replacement[i] = renumber_atoms[int(replacement[i])]
        out.append('{:>6}{:>6}{:>6}{:>6}{:>13}{:>13}\n'.format(*tuple(replacement[:6])))
        angles_section += 1

    ################# DIHEDRALS #####################

    out.append('\n[ dihedrals ]\n')
    dihedrals_section = 0

    while a[dihedrals_section].count('[ dihedrals ]') == 0:
        dihedrals_section += 1

    dihedrals_section += 2
    out.append(a[dihedrals_section])
    dihedrals_section += 1

    while a[dihedrals_section] != '\n':
        replacement = a[dihedrals_section].split()
        for i in range(4):
            replacement[i] = renumber_atoms[int(replacement[i])]
        out.append('{:>6}{:>6}{:>6}{:>6}{:>6}{:>11}{:>11}{:>11}{:>11}{:>11}{:>11}\n'.format(*tuple(replacement[:11])))
        dihedrals_section += 1

    ################ IMPROPER DIHEDRALS ###############

    out.append('\n[ dihedrals ] ; impropers\n')
    dihedrals_imp_section = 0

    while a[dihedrals_imp_section].count('[ dihedrals ] ; impropers') == 0:
        dihedrals_imp_section += 1

    dihedrals_imp_section += 2
    out.append(a[dihedrals_imp_section])
    dihedrals_imp_section += 1

    while a[dihedrals_imp_section] != '\n':
        replacement = a[dihedrals_imp_section].split()
        for i in range(4):
            replacement[i] = renumber_atoms[int(replacement[i])]
        out.append('{:>6}{:>6}{:>6}{:>6}{:>6}{:>9}{:>10}{:>4}\n'.format(*tuple(replacement[:8])))
        dihedrals_imp_section += 1

    # The rest can stay
    remainder = 0
    while a[remainder].count('[ system ]') == 0:
        remainder += 1

    for i in range(remainder - 1, len(a)):
        out.append(a[i])

    with open(args.output_top, 'w') as f:

        for line in out:

            f.write(line)

    # now rewrite .gro file

    a = []
    with open(args.gro, 'r') as f:

        for line in f:
            a.append(line)

    out = []
    out.append(a[0])
    out.append(a[1])

    for i in range(1, atom_no):
        k = list(renumber_atoms.keys())[list(renumber_atoms.values()).index(i)]
        replacement = a[k + 1].split()  # add two for top lines
        replacement[0] = renumber_res[k]  # residue number associated with that atom number
        replacement[1] = used_residues[renumber_res[k] - 1]
        replacement[3] = i

        out.append('{:>5}{:>5}{:>5}{:>5}{:>8}{:>8}{:>8}\n'.format(*tuple(replacement)))

    out.append(a[-1])

    with open(args.output_gro, 'w') as f:

        for line in out:
            f.write(line)