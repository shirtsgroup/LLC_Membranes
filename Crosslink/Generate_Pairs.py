#!/usr/bin/python

# Script to generate a pairs list based on bond information in a Gromacs topology
# The pairs list lists atoms that are 1-4 neighbors
import numpy as np
import os
import sys
import time
import pprint as pp
from scipy.sparse import csr_matrix
import genpairs

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

bonds_list = []
for i in range(bonds_index + 2, bond_count):
    bonds_list.append(a[i])


# [ dihedrals ] ; propers

dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
while a[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
    dihedrals_p_index += 1

ndp = 0  # number of lines in the 'dihedrals ; proper' section
dihedrals_p_count = dihedrals_p_index + 3  # keep track of index of a
dihedrals_prev = []
while a[dihedrals_p_count] != '\n':
    dihedrals_prev.append([int(a[dihedrals_p_count][0:6]), int(a[dihedrals_p_count][6:13]), int(a[dihedrals_p_count][13:20]), int(a[dihedrals_p_count][20:27])])
    dihedrals_p_count += 1
    ndp += 1

# [ dihedrals ] ; impropers

dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
while a[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
    dihedrals_imp_index += 1

ndimp = 0  # number of lines in the 'dihedrals ; impropers' section
dihedrals_imp_count = dihedrals_imp_index + 2
for i in range(dihedrals_imp_count, len(a)):  # This is the last section in the input .itp file
    dihedrals_imp_count += 1
    ndimp += 1

# [ pairs ]

pairs_index = 0  # find index where [ pairs ] section begins
while a[pairs_index].count('[ pairs ]') == 0:
    pairs_index += 1

npair = 0  # number of lines in the 'pairs' section
pairs_count = pairs_index + 2  # keep track of index of a
pairs_list_prev = []
while a[pairs_count] != '\n':
    pairs_list_prev.append([int(a[pairs_count][0:6]), int(a[pairs_count][6:13])])  # record pairs that already exist
    pairs_count += 1
    npair += 1

# [ angles ]

angles_index = 0  # find index where [ angles ] section begins
while a[angles_index].count('[ angles ]') == 0:
    angles_index += 1

na = 0  # number of lines in the 'angles' section
angle_count = angles_index + 2  # keep track of index of a
angles_prev = []
while a[angle_count] != '\n':
    angles_prev.append([int(a[angle_count][0:6]), int(a[angle_count][6:13]), int(a[angle_count][13:20])])
    angle_count += 1
    na += 1
print na
start = time.time()

atoms_of_interest = []
for i in range(0, 137):
    atoms_of_interest.append(i + 1)


# Make a condensed list of bonds including only the bonds which we are interested in as well as a full list of bonds:

def read_lists(a, atoms_of_interest, start, end):
    list = []
    condensed = []
    for i in range(start, end):  # extract atom numbers in [ bonds ] section in order of appearance
        list.append([int(a[i][0:6]), int(a[i][6:13])])
        if int(a[i][0:6]) in atoms_of_interest:
            condensed.append([int(a[i][0:6]), int(a[i][6:13])])
        if int(a[i][6:13]) in atoms_of_interest:
            condensed.append([int(a[i][0:6]), int(a[i][6:13])])

    uniq_condensed = []
    for i in range(0, len(condensed)):
        a = condensed[i]
        b = [condensed[i][1], condensed[i][0]]
        if a not in uniq_condensed:
            if b not in uniq_condensed:
                uniq_condensed.append(a)

    return list, uniq_condensed

list, uniq_condensed = read_lists(a, atoms_of_interest, bonds_index + 2, bond_count)

start = time.time()


def gen_pairs_list(list, condensed):
    pairs = []
    angles = []
    dihedrals = []
    for i in range(0, len(condensed)):  # look at all bonds
        neighbor3 = []  # redefine lists for each loop since each loops looks at a different bond
        neighbor4 = []
        for j in range(0, len(list)):  # compare to all other bonds
            if list[j][0] == condensed[i][1]:  # look for 1-3 pairs by comparing first entry of each row of [ pairs ] section
                neighbor3.append(list[j][1])
                angles.append([condensed[i][0], list[j][0], list[j][1]])
            elif list[j][1] == condensed[i][1]:  # look for 1-3 pairs again
                if list[j][0] != condensed[i][0]:
                    neighbor3.append(list[j][0])
                    angles.append([condensed[i][0], list[j][1], list[j][0]])
        for h in neighbor3:
            for l in range(0, len(list)):
                if list[l][0] == h:
                    if list[l][1] != condensed[i][1]:
                        dihedrals.append([condensed[i][0], condensed[i][1], h, list[l][1]])
                        neighbor4.append(list[l][1])
                elif list[l][1] == h:
                    if list[l][0] != condensed[i][1]:
                        dihedrals.append([condensed[i][0], condensed[i][1], h, list[l][0]])
                        neighbor4.append(list[l][0])
        count = 0
        for m in neighbor4:
            pairs.append([condensed[i][0], neighbor4[count]])
            count += 1

        neighbor3 = []
        neighbor4 = []
        for j in range(0, len(list)):
            if list[j][1] == condensed[i][0]:
                angles.append([condensed[i][1], list[j][1], list[j][0]])
                neighbor3.append(list[j][0])
            elif list[j][0] == condensed[i][0]:
                if list[j][1] != condensed[i][1]:
                    angles.append([condensed[i][1], list[j][0], list[j][1]])
                    neighbor3.append(list[j][1])
        for h in neighbor3:
            for l in range(0, len(list)):
                if list[l][0] == h:
                    if list[l][1] != condensed[i][0]:
                        dihedrals.append([condensed[i][1], condensed[i][0], h, list[l][1]])
                        neighbor4.append(list[l][1])
                elif list[l][1] == h:
                    if list[l][0] != condensed[i][0]:
                        dihedrals.append([condensed[i][1], condensed[i][0], h, list[l][0]])
                        neighbor4.append(list[l][0])

        count = 0
        for m in neighbor4:
            pairs.append([condensed[i][1], neighbor4[count]])
            count += 1
    return pairs, dihedrals, angles


end1 = time.time()
print end1 - start

def uniq_pair(pairs, pairs_list_prev):
    for i in range(0, len(pairs)):
        a = pairs[i]  # the order it is read from the pairs list
        b = [pairs[i][1], pairs[i][0]]
        if a not in pairs_list_prev:
            if b not in pairs_list_prev:
                pairs_list_prev.append(a)
    return pairs_list_prev

def uniq_dihedral(dihedrals, dihedrals_prev):
    for i in range(0, len(dihedrals)):
        a = dihedrals[i]
        b = [dihedrals[i][3], dihedrals[i][2], dihedrals[i][1], dihedrals[i][0]]
        if a not in dihedrals_prev:
            if b not in dihedrals_prev:
                dihedrals_prev.append(dihedrals[i])
    return dihedrals_prev

def uniq_angles(angles, angles_prev):
    for i in range(0, len(angles)):
        a = angles[i]
        b = [angles[i][2], angles[i][1], angles[i][0]]
        if a not in angles_prev:
            if b not in angles_prev:
                angles_prev.append(angles[i])
    return angles_prev

pairs, dihedrals, angles = gen_pairs_list(list, uniq_condensed)

unique_pairs = uniq_pair(pairs, pairs_list_prev)
unique_dihedrals = uniq_dihedral(dihedrals, dihedrals_prev)
unique_angles = uniq_angles(angles, angles_prev)
print len(unique_pairs)
print len(unique_dihedrals)
print len(unique_angles)
# for i in range(0, len(unique_pairs)):
#     print '{:4d}{:4d}'.format(unique_pairs[i][0], unique_pairs[i][1])
#
# print len(unique_pairs)
#
# unique_pairs = uniq_pair(pairs1, pairs2, pairs_list_prev)

