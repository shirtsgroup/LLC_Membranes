#!/usr/bin/python
# Script to make a .itp for an entire monomer assembly

import numpy as np
import math
import sys
import os
import argparse

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  #Directory where this file is located


def initialize():
    parser = argparse.ArgumentParser(description = 'Build LLC Structure')

    parser.add_argument('-t', '--type', default='LLC', help = 'Type of monomer being used')
    parser.add_argument('-l', '--layers', default=20, help = 'Number of layers')
    parser.add_argument('-P', '--no_pores', default=4, help = 'Number of Pores')
    parser.add_argument('-o', '--no_monomers', default=6, help = 'Number of monomers per layer')
    parser.add_argument('-x', '--xlink', default='on', help = 'Whether this is going to be cross linked')
    parser.add_argument('-m', '--monomer', default='monomer2', help = 'monomer being used')
    parser.add_argument('-O', '--output', default='assembly.itp', help= 'Name of output file')

    args = parser.parse_args()
    return args


def write_file(a, xlink, output, no_mon):
    # print up to ' [ atoms ] ' since everything before it does not need to be modified

    atoms_index, bonds_index, pairs_index, angles_index, dihedrals_p_index, \
    dihedrals_imp_index, vsite_index = get_indices(a, xlink)

    f = open('%s' % output, 'w')

    for i in range(0, atoms_index + 2):  # prints up to and including [ atoms ] in addition to the header line after it
        f.write(a[i])

    # [ atoms ]

    atoms_count = atoms_index + 2
    nr = 0  # number of atoms
    while a[atoms_count] != '\n':
        atoms_count += 1  # increments the while loop
        nr += 1  # counts number of atoms

    for i in range(0, no_mon):  # print atom information for each monomer
        for k in range(0, nr):  # getting the number right
            f.write('{:5d}{:25s}{:5d}{:}'.format(i*nr + k + 1, a[k + atoms_index + 2][6:29],
                                               i*nr + int(a[k + atoms_index + 2][29:34]),
                                               a[k + atoms_index + 2][34:len(a[k + atoms_index + 2])]))

    f.write("\n")  # space in between sections

    # [ bonds ]

    f.write(a[bonds_index] + a[bonds_index + 1])

    nb = 0  # number of lines in the 'bonds' section
    bond_count = bonds_index + 2
    while a[bond_count] != '\n':
        bond_count += 1  # increments while loop
        nb += 1  # counting number of lines in 'bonds' section

    for i in range(0, no_mon):
        for k in range(0, nb):
            f.write('{:6d}{:7d}{:}'.format(i*nr + int(a[k + bonds_index + 2][0:6]), i*nr + int(a[k + bonds_index + 2][6:14]),
                                         a[k + bonds_index + 2][14:len(a[k+ atoms_index + 2])]))

    f.write("\n")  # space in between sections

    # [ pairs ]

    f.write(a[pairs_index] + a[pairs_index + 1])

    npair = 0  # number of lines in the 'pairs' section
    pairs_count = pairs_index + 2  # keep track of index of a
    while a[pairs_count] != '\n':
        pairs_count += 1
        npair += 1

    for i in range(0, no_mon):
        for k in range(0, npair):
            f.write('{:6d}{:7d}{:}'.format(i*nr + int(a[k + pairs_index + 2][0:6]), i*nr + int(a[k + pairs_index + 2][6:14]),
                                         a[k + pairs_index + 2][14:len(a[k + pairs_index + 2])]))

    f.write("\n")  # space in between sections

    # [ angles ]

    f.write(a[angles_index] + a[angles_index + 1])

    na = 0  # number of lines in the 'angles' section
    angle_count = angles_index + 2  # keep track of index of a
    while a[angle_count] != '\n':
        angle_count += 1
        na += 1

    for i in range(0, no_mon):
        for k in range(0, na):
            f.write('{:6d}{:7d}{:7d}{:}'.format(i*nr + int(a[k + angles_index + 2][0:6]), i*nr + int(a[k + angles_index + 2][6:14]),
                                              i*nr + int(a[k + angles_index + 2][14:22]),
                                                         a[k + angles_index + 2][22:len(a[k + angles_index + 2])]))

    f.write("\n")  # space in between sections

    # [ dihedrals ] ; propers

    f.write(a[dihedrals_p_index] + a[dihedrals_p_index + 2])  # +2 because there is extra info that we don't need on one line

    ndp = 0  # number of lines in the 'dihedrals ; proper' section
    dihedrals_p_count = dihedrals_p_index + 3  # keep track of index of a
    while a[dihedrals_p_count] != '\n':
        dihedrals_p_count += 1
        ndp += 1

    for i in range(0, no_mon):
        for k in range(0, ndp):
            f.write('{:6d}{:7d}{:7d}{:7d}{:}'.format(i*nr + int(a[k + dihedrals_p_index + 3][0:6]),
                                                   i*nr + int(a[k + dihedrals_p_index + 3][6:14]),
                                                   i*nr + int(a[k + dihedrals_p_index + 3][14:22]),
                                                   i*nr + int(a[k + dihedrals_p_index + 3][22:30]),
                                                   a[k + dihedrals_p_index + 3][30:len(a[k + dihedrals_p_index + 3])]))

    f.write("\n")  # space in between sections

    # [ dihedrals ] ; impropers

    f.write(a[dihedrals_imp_index] + a[dihedrals_imp_index + 2]),
    ndimp = 0  # number of lines in the 'dihedrals ; impropers' section
    dihedrals_imp_count = dihedrals_imp_index + 3

    while a[dihedrals_imp_count] != '\n':
        dihedrals_imp_count += 1
        ndimp += 1

    # Can't have any space at the bottom of the file for this loop to work
    for i in range(0, no_mon):
        for k in range(0, ndimp):
            f.write('{:6d}{:7d}{:7d}{:7d}{:}'.format(i*nr + int(a[k + dihedrals_imp_index + 3][0:6]),
                                                   i*nr + int(a[k + dihedrals_imp_index + 3][6:14]),
                                                   i*nr + int(a[k + dihedrals_imp_index + 3][14:22]),
                                                   i*nr + int(a[k + dihedrals_imp_index + 3][22:30]),
                                                   a[k + dihedrals_imp_index + 3][30:len(a[k + dihedrals_imp_index + 3])]))
    f.write("\n")  # space in between sections

    # [ virtual_sites4 ]
    if xlink == 'on':
        f.write(a[vsite_index] + a[vsite_index + 1]),
        nv = 0
        vsite_count = vsite_index + 2

        for i in range(vsite_count, len(a)):  # This is the last section in the input .itp file
            vsite_count += 1
            nv += 1

        # Make sure there is no space at the bottom of the topology if you are getting errors
        for i in range(0, no_mon):
            for k in range(0, nv):
                f.write('{:<8d}{:<6d}{:<6d}{:<6d}{:<8d}{:<8d}{:<11}{:<11}{:}'.format(i*nr + int(a[k + vsite_index + 2][0:8]),
                                                       i*nr + int(a[k + vsite_index + 2][8:14]),
                                                       i*nr + int(a[k + vsite_index + 2][14:20]),
                                                       i*nr + int(a[k + vsite_index + 2][20:26]),
                                                       i*nr + int(a[k + vsite_index + 2][26:34]),
                                                       int(a[k + vsite_index + 2][34:42]), a[k + vsite_index + 2][42:53],
                                                       a[k + vsite_index + 2][53:64],
                                                       a[k + vsite_index + 2][64:len(a[k + vsite_index + 2])]))
    f.close()


def get_indices(a, xlink):
    # find the indices of all fields that need to be modified
    atoms_index = 0  # find index where [ atoms ] section begins
    while a[atoms_index].count('[ atoms ]') == 0:
        atoms_index += 1

    bonds_index = 0  # find index where [ bonds ] section begins
    while a[bonds_index].count('[ bonds ]') == 0:
        bonds_index += 1

    pairs_index = 0  # find index where [ pairs ] section begins
    while a[pairs_index].count('[ pairs ]') == 0:
        pairs_index += 1

    angles_index = 0  # find index where [ angles ] section begins
    while a[angles_index].count('[ angles ]') == 0:
        angles_index += 1

    dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
    while a[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
        dihedrals_p_index += 1

    dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
    while a[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
        dihedrals_imp_index += 1

    if xlink == 'on':
        vsite_index = 0  # find index where [ dihedrals ] section begins (propers)
        while a[vsite_index].count('[ virtual_sites4 ]') == 0:
            vsite_index += 1
    else:
        vsite_index = 0

    return atoms_index, bonds_index, pairs_index, angles_index, dihedrals_p_index, dihedrals_imp_index, vsite_index

if __name__ == "__main__":
    args = initialize()
    no_mon = args.layers * args.no_pores * args.no_monomers  # How many monomer will be in this .itp file

    if args.type == 'LLC':
        if args.xlink == 'on':
            itp = '%s_dummy.itp' % args.monomer
        else:
            itp = '%s.itp' % args.monomer
    elif args.type == 'BCC':
        itp = 'BCC_mon.itp'

    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory where this file is located
    f = open("%s/../Structure-Files/Monomer_Tops/%s" % (location, itp), "r")

    a = []
    for line in f:
        a.append(line)

    f.close()

    write_file(a, '%s' % args.xlink, '%s' % args.output)



