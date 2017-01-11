#!/usr/bin/python

"""
    Converts a gro file to a .xyz file
"""

import argparse
import LC_class
import numpy as np


def initialize():
    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-f', '--file', default='wiggle.gro', help='Name of file to read')
    parser.add_argument('-o', '--output', default='wiggle.xyz', help='Name of file to be output')
    parser.add_argument('-l', '--lc', default='HII', help='Type of liquid crystal')
    args = parser.parse_args()
    return args


def coords_and_identity(file, lc):

    exec 'residues = LC_class.%s.residues' % lc

    # find the first line containing coordinates
    coord_start = 0  # start at the top
    while str.strip(file[coord_start][5:10]) not in residues:
        coord_start += 1

    no_atoms = len(file) - coord_start - 1  # subtract 1 due to the last line of the file containing box vectors

    coords = np.zeros([3, no_atoms])
    identity = np.zeros([no_atoms], dtype=object)

    for i in range(coord_start, no_atoms + coord_start):
        index = i - coord_start
        coords[0, index] = float(file[i][21:29])  # Use this to read specific entries in a text file
        coords[1, index] = float(file[i][29:38])  # makes sure I backtrack far enough to get all digits(i.e.38 instead of 42)
        coords[2, index] = float(file[i][38:45])
        identity[index] = str(file[i][9:15])  # hold name of atom (C, H or O)
        if identity[index].count('C') == 1:
            identity[index] = 'C'
        elif identity[index].count('H') == 1:
            identity[index] = 'H'
        elif identity[index].count('O') == 1:
            identity[index] = 'O'
        elif identity[index].count('NA') == 1:
            identity[index] = 'NA'
        else:
            print 'unrecognized atoms(s)'

    return coords, identity


def write_xyz(filename, identity, coords):

    f = open('%s' % filename, 'w')

    for i in range(len(identity)):
        f.write('{:>2}{:12}{:12}{:12}'.format(identity[i], coords[0, i], coords[1, i], coords[2, i]) + '\n')

    f.close()

if __name__ == "__main__":
    args = initialize()

    f = open('%s' % args.file, 'r')  # .gro file whose atomic coordinates will be read
    a = []  # list to hold lines of file
    for line in f:
        a.append(line)
    f.close()

    coords, identity = coords_and_identity(a, '%s' % args.lc)

    write_xyz('%s' % args.output, identity, coords)