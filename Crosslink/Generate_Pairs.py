#!/usr/bin/python

# Script to generate a pairs list based on bond information in a Gromacs topology
# The pairs list lists atoms that are 1-4 neighbors

import numpy as np
import os
import sys

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

# locate [ bonds ] section

bonds_index = 0  # find index where [ bonds ] section begins
while a[bonds_index].count('[ bonds ]') == 0:
    bonds_index += 1

nb = 0  # number of lines in the 'bonds' section
bond_count = bonds_index + 2
while a[bond_count] != '\n':
    bond_count += 1  # increments while loop
    nb += 1  # counting number of lines in 'bonds' section

# Start finding 1-4 pairs

pairs = []
for i in range(bonds_index + 2, bond_count):  # look at each bond in [ bonds ] section
    for k in range(1, nr + 1):  # look for all atoms
        neighbor3 = []  # 1-3 neighbors
        neighbor4 = []  # 1-4 neighbors
        if int(a[i][0:6]) == k:
            neighbor1 = int(a[i][6:13])  # This the atom's immediate neighbor which it is bonded to
            for j in range(bonds_index + 2, bond_count):
                if int(a[j][0:6]) == neighbor1:
                    if int(a[j][6:13]) != k:
                        neighbor3.append(int(a[j][6:13]))
                elif int(a[j][6:13]) == neighbor1:
                    if int(a[j][0:6]) != k:
                        neighbor3.append(int(a[j][0:6]))
            for h in neighbor3:
                for l in range(bonds_index + 2, bond_count):
                    if int(a[l][0:6]) == h:
                        if int(a[l][6:13]) != neighbor1:
                            neighbor4.append(int(a[l][6:13]))
                    elif int(a[l][6:13]) == h:
                        if int(a[l][0:6]) != neighbor1:
                            neighbor4.append(int(a[l][0:6]))
        if int(a[i][6:13]) == k:
            neighbor1 = int(a[i][0:6])  # This the atom's immediate neighbor which it is bonded to
            for j in range(bonds_index + 2, bond_count):
                if int(a[j][0:6]) == neighbor1:
                    if int(a[j][6:13]) != k:
                        neighbor3.append(int(a[j][6:13]))
                elif int(a[j][6:13]) == neighbor1:
                    if int(a[j][0:6]) != k:
                        neighbor3.append(int(a[j][0:6]))
            for h in neighbor3:
                for l in range(bonds_index + 2, bond_count):
                    if int(a[l][0:6]) == h:
                        if int(a[l][6:13]) != neighbor1:
                            neighbor4.append(int(a[l][6:13]))
                    elif int(a[l][6:13]) == h:
                        if int(a[l][0:6]) != neighbor1:
                            neighbor4.append(int(a[l][0:6]))
            # print neighbor4
        count = 0
        for m in neighbor4:
            pairs.append('{:6d}{:7d}{:7d}'.format(k, neighbor4[count], 1))
            count += 1

# eliminate duplicates

def uniq_pair(pairs):
    output = []
    for i in range(0, len(pairs)):
        a = pairs[i]  # the order it is read from the pairs list
        b = (len(pairs[i][6:13]) - len(str.strip(pairs[i][6:13])) - 1)*' ' + str.strip(pairs[i][6:13]) + ' ' + \
            pairs[i][0:6] + pairs[i][13:len(pairs[i])]  # switching the order of the pairs for comparison
        if a not in output:
            if b not in output:
                output.append(pairs[i])
    return output

unique_pairs = uniq_pair(pairs)
# Now print them in order

ordered_unique_pairs = []
for i in range(1, nr + 1):
    no_same_atom = []
    for k in range(0, len(unique_pairs)):
        if int(unique_pairs[k][0:6]) == i:
            no_same_atom.append(unique_pairs[k])
    for l in range(1, nr + 1):
        for j in range(0, len(no_same_atom)):
            if int(no_same_atom[j][6:13]) == l:
                ordered_unique_pairs.append(no_same_atom[j])

# for i in range(0, len(ordered_unique_pairs)):
#     print ordered_unique_pairs[i]

# TESTING -- it works! Boom!

pairs_index = 0  # find index where [ pairs ] section begins
while a[pairs_index].count('[ pairs ]') == 0:
    pairs_index += 1

npair = 0  # number of lines in the 'pairs' section
pairs_count = pairs_index + 2  # keep track of index of a
while a[pairs_count] != '\n':
    pairs_count += 1
    npair += 1

print a[pairs_index + 2][0:20]
print a[pairs_count][0:20]

test = []
for i in range(pairs_index + 2, pairs_count):
    test.append(a[i][0:20])

count = 0
for i in range(0, len(test)):
    a = test[i]  # the order it is read from the pairs list
    b = (len(test[i][6:13]) - len(str.strip(test[i][6:13])) - 1)*' ' + str.strip(test[i][6:13]) + ' ' + \
        test[i][0:6] + test[i][13:len(test[i])]  # switching the order of the pairs for comparison
    if a or b in ordered_unique_pairs:
        count += 1

print count
