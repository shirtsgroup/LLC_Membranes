#!/usr/bin/env python

"""
Apply flat-bottomed position restraints based on a reference position
"""

import numpy as np
import argparse
import mdtraj as md
from scipy import spatial
import os


def initialize():

    parser = argparse.ArgumentParser(description='Apply flat-bottomed position restraints based on a reference position')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='A coordinate file')
    parser.add_argument('-c', '--coordinate', default=[4, 4, 4], help='Reference coordinate for position restraint')
    parser.add_argument('--exclude', action="store_true", help='Exclude things from region rather than confine them')
    parser.add_argument('-r', '--radius', default=1, help='radius or distance from reference coordinate defining the region')
    parser.add_argument('-f', '--function', default=2, help='type of region. 1 = sphere, 2 = cylinder etc. (see gromacs docs)')
    parser.add_argument('-o', '--output', default='restraints.itp', help='name of output file')
    parser.add_argument('-k', '--force', default=10, help='Force constant')

    args = parser.parse_args()

    return args


def get_indices(a):
    # find the indices of all fields that need to be modified
    atoms_index = 0  # find index where [ atoms ] section begins
    while a[atoms_index].count('[ atoms ]') == 0:
        atoms_index += 1

    bonds_index = 0  # find index where [ bonds ] section begins
    while a[bonds_index].count('[ bonds ]') == 0:
        bonds_index += 1

    exclusions_index = 0  # find index where [ pairs ] section begins
    while a[exclusions_index].count('[ exclusions ]') == 0:
        exclusions_index += 1

    angles_index = 0  # find index where [ angles ] section begins
    while a[angles_index].count('[ angles ]') == 0:
        angles_index += 1

    settles_index = 0  # find index where [ angles ] section begins
    while a[settles_index].count('[ settles ]') == 0:
        settles_index += 1

    return atoms_index, bonds_index, angles_index, settles_index, exclusions_index


def write_assembly(itp, output, n):

    """
    :param itp: Name of .itp file for single component
    :param n: number of molecules in the assembly
    :param output: name of output file
    :return:
    """

    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Location of this script

    with open("%s/../top/Forcefields/gaff/%s" % (location, 'tip3p.itp'), "r") as f:

        a = []
        for line in f:
            a.append(line)

    atoms_index, bonds_index, angles_index, settles_index, exclusions_index = get_indices(a)
    f = open('%s' % output, 'w')

    # print up to ' [ atoms ] ' since everything before it does not need to be modified
    for i in range(0, atoms_index + 2):  # prints up to and including [ atoms ] in addition to the header line after it
        f.write(a[i])

    # [ atoms ]

    atoms_count = atoms_index + 2
    nr = 0  # number of atoms
    while a[atoms_count] != '\n':
        atoms_count += 1  # increments the while loop
        nr += 1  # counts number of atoms

    for i in range(0, n):  # print atom information for each monomer
        for k in range(0, nr):  # getting the number right
            f.write('{:6d}{:}{:5s}{:5d}{:}'.format(i*nr + k + 1, ' ', a[k + atoms_index + 2][6:18],
                                               i*nr + k + 1, a[k + atoms_index + 2][24:len(a[k + atoms_index + 2])]))

    f.write("\n#ifndef FLEXIBLE\n\n")

    # [ settles ]

    f.write(a[settles_index] + a[settles_index + 1])

    ns = 0  # number of lines in the 'bonds' section
    settles_count = settles_index + 2
    while a[settles_count] != '\n':
        settles_count += 1  # increments while loop
        ns += 1  # counting number of lines in 'bonds' section

    for i in range(0, n):
        for k in range(0, ns):
            f.write('{:6d}{:6d}{:}'.format(i*nr + int(a[k + settles_index + 2][:7]), i*nr + int(a[k + settles_index + 2][7:15]),
                                         a[k + settles_index + 2][15:len(a[k + settles_index + 2])]))

    f.write('\n')

    # [ exclusions ]

    f.write(a[exclusions_index])

    ne = 0  # number of lines in the 'bonds' section
    exclusions_count = exclusions_index + 1
    while a[exclusions_count] != '\n':
        exclusions_count += 1  # increments while loop
        ne += 1  # counting number of lines in 'bonds' section

    for i in range(0, n):
        for k in range(0, ne):
            f.write('{:6d}{:6d}{:6d}\n'.format(i*nr + int(a[k + exclusions_index + 1][:7]), i*nr + int(a[k + exclusions_index + 1][7:15]),
                                         i*nr + int(a[k + exclusions_index + 1][15:20])))

    f.write('\n#else\n\n')

    # [ bonds ]

    f.write(a[bonds_index] + a[bonds_index + 1])

    nb = 0  # number of lines in the 'bonds' section
    bond_count = bonds_index + 2
    while a[bond_count] != '\n':
        bond_count += 1  # increments while loop
        nb += 1  # counting number of lines in 'bonds' section
    nb -= 1

    for i in range(0, n):
        for k in range(0, nb):
            f.write('{:6d}{:7d}{:}'.format(i*nr + int(a[k + bonds_index + 2][:7]), i*nr + int(a[k + bonds_index + 2][7:16]),
                                         a[k + bonds_index + 2][16:len(a[k + atoms_index + 2])]))

    f.write("\n")  # space in between sections

    # [ angles ]

    f.write(a[angles_index] + a[angles_index + 1])

    na = 0  # number of lines in the 'angles' section
    angle_count = angles_index + 2  # keep track of index of a
    while a[angle_count] != '\n':
        angle_count += 1
        na += 1

    for i in range(0, n):
        for k in range(0, na):
            f.write('{:6d}{:7d}{:7d}{:}'.format(i*nr + int(a[k + angles_index + 2][0:6]), i*nr + int(a[k + angles_index + 2][6:14]),
                                              i*nr + int(a[k + angles_index + 2][14:22]),
                                                         a[k + angles_index + 2][22:len(a[k + angles_index + 2])]))

    f.write("\n#endif\n")  # space in between sections

    f.close()




if __name__ == "__main__":

    args = initialize()

    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Location of this script

    g = args.function
    if args.exclude:
        r = -args.radius
    else:
        r = args.radius
    k = args.force

    t = md.load(args.gro)

    pos = t.xyz[0, :, :]

    nn = spatial.KDTree(pos).query(args.coordinate)[1]  # atom closest to coordinate

    write_assembly('tip3p.itp', args.output, t.n_residues)

    with open(args.output, 'a') as f:

        f.write('\n[ position_restraints ]\n')

        for i in range(pos.shape[0]):
            if i != nn:
                f.write('{:6d}{:6d}{:1s}{:3d}{:1s}{:3d}{:1s}{:5f}\n'.format(i + 1, nn,'', g, '', r, '', k))