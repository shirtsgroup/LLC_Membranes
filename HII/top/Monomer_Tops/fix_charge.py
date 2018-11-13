#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
from past.utils import old_div
import argparse
import numpy as np


def initialize():

    parser = argparse.ArgumentParser(description='Normalize charge so that it sums to the system formal charge')

    parser.add_argument('-t', '--top', type=str, help='name of topology file '
                                                      '(.itp)')
    parser.add_argument('-c', '--charge', default=-1, type=int, help='formal charge on molecule')
    parser.add_argument('-o', '--output', default='out.top', type=str, help='Name of output file')

    args = parser.parse_args()

    return args


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

if __name__ == "__main__":

    args = initialize()

    f = open(args.top, 'r')

    a = []
    for line in f:
        a.append(line)

    # find the [ atoms ] section
    start = 0
    while a[start].count('[ atoms ]') == 0:
        start += 1

    count = start + 2
    charges = []
    while a[count] != '\n':
        charges.append(float(a[count].split()[6]))
        count += 1

    while abs(sum(charges) - args.charge) > 10**-9:
        tot = sum(charges)
        if args.charge != 0:
            charges = [round(old_div(i,abs(tot)), 6) for i in charges]
        else:
            mean = np.mean(charges)
            charges = [i - mean for i in charges]

    new_charges = replacecharge(a, charges)

    f = open(args.output, 'w')
    for line in new_charges:
        f.write(line)

    f.close()
