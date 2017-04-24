#!/usr/bin/python

"""
Use this to take the charges out of a mol2 file and assign them in place of
existing charges in a .itp file
"""

import argparse
import numpy as np


def initialize():

    parser = argparse.ArgumentParser(description='Replace charges in .itp file with charges in mol2 file')

    parser.add_argument('-m', '--mol2', default='out.mol2', help='Name of mol2 file')
    parser.add_argument('-t', '--itp', default='monomer2.itp', help='Name of topology file in which charges will be replaced')
    parser.add_argument('-o', '--out', default='new_charges.out', help='Name of output file')
    parser.add_argument('-c', '--charge', default=-1, type=int, help='The formal charge of the whole molecule')

    args = parser.parse_args()

    return args


def extractcharge(mol2, formalcharge):
    """
    :param mol2: the input mol2 file in list format
    """
    # find the line where the atoms start being listed
    atoms = 0
    while mol2[atoms].count('ATOM') == 0:
        atoms += 1

    charges = []
    atoms += 1
    while mol2[atoms].count('@') == 0:

        charges.append(float(mol2[atoms][-10:(len(mol2[atoms]) - 1)]))
        atoms += 1

    charge = np.array(charges)

    # the total charge on the system needs to be corrected to the normalized charge

    correction = (float(formalcharge) - sum(charge)) / len(charge)

    for i in range(len(charge)):
        charge[i] += correction

    return charge


def replacecharge(top, newcharges):
    """
    :param top: the input topology file in list format
    """
    # find the line where the [ atoms ] section begins
    atoms = 0
    while top[atoms].count('[ atoms ]') == 0:
        atoms += 1

    atoms += 1
    count = 0
    while count < len(newcharges):
        if top[atoms][0] == ";":
            atoms += 1
        else:
            a = top[atoms].split()
            a[6] = newcharges[count]
            top[atoms] = "{:6d}{:>5s}{:6d}{:>6s}{:>6s}{:5d}{:13f}{:13f}\n".format(int(a[0]), a[1], int(a[2]), a[3], a[4], int(a[5]), float(a[6]), float(a[7]))
            count += 1
            atoms += 1

    return top


if __name__ == '__main__':

    args = initialize()

    # read input files and put each line into a list
    f = open('%s' % args.mol2, 'r')
    mol2 = []
    for line in f:
        mol2.append(line)
    f.close()

    f = open('%s' % args.itp, 'r')
    top = []
    for line in f:
        top.append(line)
    f.close()

    charges = extractcharge(mol2, args.charge)

    new_charges = replacecharge(top, charges)

    f = open('%s' % args.out, 'w')
    for line in new_charges:
        f.write(line)
    f.close()
