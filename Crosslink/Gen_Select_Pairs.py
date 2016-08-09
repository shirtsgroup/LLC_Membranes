#!/usr/bin/python

# Script to generate a pairs list based on bond information in a Gromacs topology
# The pairs list lists atoms that are 1-4 neighbors

import numpy as np
import os
import sys
import time
import pprint as pp
from scipy.sparse import csr_matrix

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

f = open('%s/../Structure-Files/Monomer_Tops/HII_mon.itp' %location, 'r')

a = []
for line in f:
    a.append(line)

# locate [ atoms ] section

atoms_index = 0  # find index where [ atoms ] section begins
while a[atoms_index].count('[ atoms ]') == 0:
    atoms_index += 1

# count number of atoms
atoms_count = atoms_index + 2
nr = 0  # number of lines in 'atoms' section
while a[atoms_count] != '\n':
    atoms_count += 1  # increments the while loop
    nr += 1  # counts number of atoms

# Extract atom name information
atom_names = []
for i in range(atoms_index + 2, atoms_count):
    atom_names.append(str.strip(a[i][23:30]))

# Create a dictionary relating atom name to atom number

atoms_dict = dict()
for i in range(0, len(atom_names)):
    if atom_names[i] in atoms_dict:
        atoms_dict[atom_names[i]].append(i + 1)
    else:
        atoms_dict[atom_names[i]] = [i + 1]
atoms_dict['O5'].append(190)

# Define atoms of interest that belong to each tail

tail1_atoms = ['O5', 'O10', 'C18', 'C19', 'C20', 'H24', 'H25', 'H26']
tail2_atoms = ['O6', 'O9', 'C33', 'C34', 'C32', 'H49', 'H50', 'H51']
tail3_atoms = ['O7', 'O8', 'C47', 'C48', 'C46', 'H74', 'H75', 'H76']

tail1 = dict()
tail2 = dict()
tail3 = dict()

for key in atoms_dict:
    if key in tail1_atoms:
        tail1[key] = atoms_dict[key]
    if key in tail2_atoms:
        tail2[key] = atoms_dict[key]
    if key in tail3_atoms:
        tail3[key] = atoms_dict[key]


# locate [ bonds ] section

bonds_index = 0  # find index where [ bonds ] section begins
while a[bonds_index].count('[ bonds ]') == 0:
    bonds_index += 1

nb = 0  # number of lines in the 'bonds' section
bond_count = bonds_index + 2
while a[bond_count] != '\n':
    bond_count += 1  # increments while loop
    nb += 1  # counting number of lines in 'bonds' section

pairs_index = 0  # find index where [ pairs ] section begins
while a[pairs_index].count('[ pairs ]') == 0:
    pairs_index += 1

npair = 0  # number of lines in the 'pairs' section
pairs_count = pairs_index + 2  # keep track of index of a
while a[pairs_count] != '\n':
    pairs_count += 1
    npair += 1

# Make list of current 1-4 pairs
pairs_current = []
for i in range(pairs_index + 2, pairs_count):
    pairs_current.append(a[i][0:20])

# list of atoms involved in cross-linking (even if they are just 1-4 pairs with a cross-linking atom)

atoms_list = ['O5', 'O6', 'O7', '08', 'O9', 'O10', 'C32', 'C18', 'C46', 'C33', 'C19', 'C47', 'C34', 'C20', 'C47', 'H49',
              'H24', 'H74', 'H50', 'H25', 'H75', 'H51', 'H26', 'H76']

n = 27
all_numbers = [n]
if atom_names[n - 1] in tail1:
    all_numbers.append(n - 3)
    all_numbers.append(n - 2)
    all_numbers.append(n - 1)
    all_numbers.append(n + 33)
    all_numbers.append(n + 58)
    all_numbers.append(n + 59)
    all_numbers.append(n + 60)

print all_numbers
print tail1
# Given atom numbers, identify which atoms they are
atom_no = 190
for i in range(0, len(tail1.values())):
    if atom_no in tail1.values()[i]:
        print 'hello'


# TESTING
#
#
# test = []
# for i in range(pairs_index + 2, pairs_count):
#     test.append(a[i][0:20])
#
# count = 0
# for i in range(0, len(test)):
#     a = test[i]  # the order it is read from the pairs list
#     b = (len(test[i][6:13]) - len(str.strip(test[i][6:13])) - 1)*' ' + str.strip(test[i][6:13]) + ' ' + \
#         test[i][0:6] + test[i][13:len(test[i])]  # switching the order of the pairs for comparison
#     if a or b in ordered_unique_pairs:
#         count += 1
#
# if count == 349:
#     print 'It works!'