#! /bin/bash/env python

from __future__ import print_function
import argparse
import os


def initialize():

    parser = argparse.ArgumentParser(description = 'Make a dictionary of charge values for different monomers')

    # User inputs with defaults
    parser.add_argument('-m', '--mon', default='NAcarb11V', help='Name of monomer (no file extension)')

    args = parser.parse_args()

    return args


def mk_dict(monomer):

    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

    with open('%s/../top/topologies/%s.itp' % (location, monomer), 'r') as f:
        a = []
        for line in f:
            a.append(line)

    atoms_index = 0
    while a[atoms_index].count('[ atoms ]') == 0:
        atoms_index += 1

    atoms_index += 2

    charge = {}

    while a[atoms_index] != '\n':
        charge[str.strip(a[atoms_index][23:29])] = float(a[atoms_index][35:47])
        atoms_index += 1

    # f = open('dict_entries.txt', 'w')
    #
    # while a[atoms_index] != '\n':
    #     f.write("""charge['%s'] = %s \n""" %(str.strip(a[atoms_index][23:29]), float(a[atoms_index][35:47])))
    #     atoms_index += 1
    #
    # f.close()
    return charge


if __name__ == "__main__":

    args = initialize()

    charge = mk_dict(args.mon)

    print(charge)
