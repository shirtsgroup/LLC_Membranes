#!/usr/bin/python

"""
Calculate the dipole moment of a group of atoms based on their positions and partial charges
"""

import argparse
import numpy as np
#['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O4', 'O3', 'O', 'O1', 'O2']

def initialize():

    parser = argparse.ArgumentParser(description='Calculate the dipole moment of a group of atoms based on their'
                                                 'positions and partial charges')

    parser.add_argument('-c', '--coord', default='initial.gro', help='Coordinate file')
    parser.add_argument('-t', '--top', default='NAcarb11V_dummy.itp', help='Coordinate file')
    parser.add_argument('-a', '--atoms', nargs='+', default='all', help='Name of atom(s) you want to use to calculate dipole'
                                                             'moment. Specify "all" to do them all')

    args = parser.parse_args()

    return args


def dipole_moment(positions, charges):
    """
    :param positions: xyz coordinates of atoms of interest
    :param charges: charge of each atom
    :return: vector telling the direction of the dipole moment
    """

    atoms = positions.shape[0]
    dm = np.zeros(3)
    for i in range(atoms):
        dm += charges[i]*positions[i, :]

    return dm

if __name__ == "__main__":
    args = initialize()

    gro = []
    with open('%s' % args.coord, 'r') as f:
        for line in f:
            gro.append(line)

    if args.atoms == 'all':
        pos = np.zeros([len(gro) - 3, 3])
    else:
        pos = np.zeros([len(args.atoms), 3])

    atoms = pos.shape[0]

    count = 0
    for i in range(2, len(gro) - 1):
        if args.atoms == 'all':
            pos[i - 2, :] = [float(gro[i][20:28]), float(gro[i][28:36]), float(gro[i][36:44])]
        elif str.strip(gro[i][10:15]) in args.atoms:
            pos[count, :] = [float(gro[i][20:28]), float(gro[i][28:36]), float(gro[i][36:44])]
            count += 1

    top = []
    count = 0
    with open('%s' % args.top, 'r') as f:
        for line in f:
            if line.count('[ atoms ]') == 1:
                atoms_start = count
            count += 1
            top.append(line)

    charges = np.zeros(pos.shape[0])

    for i in range(atoms_start + 2, atoms + 1 + atoms_start):
        charges[i - atoms_start - 2] = float(top[i][35:48])

    dm = dipole_moment(pos, charges)

    print dm