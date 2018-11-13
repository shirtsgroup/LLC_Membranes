#!/usr/bin/env python

"""
Repartition hydrogen masses in topologies by moving mass off bonded heavy atoms.
"""

import argparse


def initialize():

    parser = argparse.ArgumentParser(description='Edit topology of a system or residue to repartition mass')

    parser.add_argument('-t', '--top', type=str, default='NAcarb11V.itp', help='Topology to be repartitioned')
    parser.add_argument('-o', '--out', type=str, default='NAcarb11V_hmr.itp', help='Output repartitioned topology')
    parser.add_argument('-w', '--weight', type=float, default=4, help='Number by which to multiply hydrogen masses')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    with open(args.top, 'r') as f:
        a = []
        for line in f:
            a.append(line)

    atom = 0
    while a[atom].count('[ atoms ]') == 0:
        atom += 1

    name = {}  # name of each atom
    mass = {}  # mass of each atom
    H = []  # serial number of all H atoms in residue
    atom += 2  # line of comments
    atom_start = atom
    while a[atom] != '\n':
        info = a[atom].split()
        mass[int(info[0])] = float(info[7])
        name[int(info[0])] = info[4]
        if 'H' in info[4]:
            H.append(info[0])
        atom += 1

    bonds = 0
    while a[bonds].count('[ bonds ]') == 0:
        bonds += 1

    bonds_start = bonds - 1 # include the space above the section
    bonds += 2
    while a[bonds] != '\n':
        info = a[bonds].split()
        if info[0] in H:
            mass[int(info[0])] *= args.weight
            mass[int(info[1])] -= ((args.weight - 1)/args.weight)*mass[int(info[0])]
        if info[1] in H:
            mass[int(info[1])] *= args.weight
            mass[int(info[0])] -= ((args.weight - 1)/args.weight)*mass[int(info[1])]
        bonds += 1

    out = []
    for i in range(atom_start):
        out.append(a[i])

    while a[atom_start] != '\n':
        info = a[atom_start].split()
        info[7] = mass[int(info[0])]
        out.append('{:>6s}{:>5s}{:>6s}{:>6s}{:>6s}{:>5s}{:>13s}{:12.6f}\n'.format(info[0], info[1], info[2], info[3],
                                                                                info[4], info[5], info[6], info[7]))
        atom_start += 1

    for i in range(bonds_start, len(a)):
        out.append(a[i])

    with open(args.out, 'w') as f:

        for line in out:
            f.write(line)