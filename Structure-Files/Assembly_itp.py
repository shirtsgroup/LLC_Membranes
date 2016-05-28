#!/usr/bin/python
# Script to make a .itp for an entire monomer assembly

import numpy as np
import math
import sys
import os
import argparse

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  #Directory where this file is located

parser = argparse.ArgumentParser(description = 'Build LLC Structure')

parser.add_argument('-t', '--type', default='LLC', help = 'Type of monomer being used')
parser.add_argument('-l', '--layers', default=20, help = 'Number of layers')
parser.add_argument('-P', '--no_pores', default=4, help = 'Number of Pores')
parser.add_argument('-o', '--no_monomers', default=6, help = 'Number of monomers per layer')

args = parser.parse_args()

if args.type == 'LLC':
    itp = 'HII_mon.itp'
elif args.type == 'BCC':
    itp = 'BCC_mon.itp'


f = open("%s/Monomer_Tops/%s" % (location, itp), "r")

a = []
for line in f:
    a.append(line)

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

# print up to ' [ atoms ] ' since everything before it does not need to be modified
for i in range(0, atoms_index + 2):  # prints up to and including [ atoms ] in addition to the header line after it
    print a[i],

# [ atoms ]

atoms_count = atoms_index + 2
nr = 0  # number of atoms
while a[atoms_count] != '\n':
    atoms_count += 1  # increments the while loop
    nr += 1  # counts number of atoms

no_mon = args.layers * args.no_pores * args.no_monomers  # How many monomer will be in this .itp file
# for i in range(0, no_mon):  # print atom information for each monomer
#     for k in range(0, nr):  # getting the number right
#         print '{:5d}{:25s}{:5d}{:}'.format(i*nr + k + 1, a[k + atoms_index + 2][6:29],
#                                            i*nr + int(a[k + atoms_index + 2][29:34]),
#                                            a[k + atoms_index + 2][34:len(a[k + atoms_index + 2])]),

print ''  # space in between sections

# [ bonds ]

print a[bonds_index],
print a[bonds_index + 1],

nb = 0  # number of lines in the 'bonds' section
bond_count = bonds_index + 2
while a[bond_count] != '\n':
    bond_count += 1  # increments while loop
    nb += 1  # counting number of lines in 'bonds' section

for i in range(0, no_mon):
    for k in range(0, nb):
        print '{:6d}{:7d}{:}'.format(i*nr + int(a[k + bonds_index + 2][0:6]), i*nr + int(a[k + bonds_index + 2][6:14]),
                                     a[k + bonds_index + 2][14:len(a[k+ atoms_index + 2])]),