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

start = time.time()

list1 = []  # list to hold atom numbers of 1st entry in each row of [ bonds ] section
list2 = []  # list to hold atom numbers of 2nd entry in each row of [ bonds ] section
for i in range(bonds_index + 2, bond_count):  # extract atom numbers in [ bonds ] section in order of appearance
    list1.append(int(a[i][0:6]))  # converts them to integers to save time later
    list2.append(int(a[i][6:13]))

xlinked_list1 = []
xlinked_list2 = []
for i in range(bonds_index + 2 + 100, bond_count):  # extract atom numbers in [ bonds ] section in order of appearance
    xlinked_list1.append(int(a[i][0:6]))  # converts them to integers to save time later
    xlinked_list2.append(int(a[i][6:13]))

xlinked_list3 = []
xlinked_list4 = []
for i in range(bonds_index + 2, bond_count - 37):  # extract atom numbers in [ bonds ] section in order of appearance
    xlinked_list3.append(int(a[i][0:6]))  # converts them to integers to save time later
    xlinked_list4.append(int(a[i][6:13]))

start = time.time()

def gen_pairs_list(list1, list2, xlinked_list1, xlinked_list2, no_bonds):
    pairs1 = []
    pairs2 = []
    dihedrals = []
    dc = 0  # dihedrals count
    for i in range(0, len(xlinked_list1)):  # look at all bonds
        neighbor3 = []  # redefine lists for each loop since each loops looks at a different bond
        neighbor4 = []
        for j in range(0, no_bonds):  # compare to all other bonds
            if list1[j] == xlinked_list2[i]:  # look for 1-3 pairs by comparing first entry of each row of [ pairs ] section
                neighbor3.append(list2[j])
            elif list2[j] == xlinked_list2[i]:  # look for 1-3 pairs again
                if list1[j] != xlinked_list1[i]:
                    neighbor3.append(list1[j])
        for h in neighbor3:
            for l in range(0, len(list1)):
                if list1[l] == h:
                    if list2[l] != xlinked_list2[i]:
                        dihedrals.append([])
                        dihedrals[dc].append(xlinked_list1[i])
                        dihedrals[dc].append(xlinked_list2[i])
                        dihedrals[dc].append(h)
                        dihedrals[dc].append(list2[l])
                        dc += 1
                        neighbor4.append(list2[l])
                elif list2[l] == h:
                    if list1[l] != xlinked_list2[i]:
                        dihedrals.append([])
                        dihedrals[dc].append(xlinked_list1[i])
                        dihedrals[dc].append(xlinked_list2[i])
                        dihedrals[dc].append(h)
                        dihedrals[dc].append(list1[l])
                        dc += 1
                        neighbor4.append(list1[l])
        count = 0
        for m in neighbor4:
            pairs1.append(xlinked_list1[i])
            pairs2.append(neighbor4[count])
            #pairs.append('{:6d}{:7d}{:7d}'.format(list1[i], neighbor4[count], 1))
            count += 1

        neighbor3 = []
        neighbor4 = []
        for j in range(0, no_bonds):
            if list2[j] == xlinked_list1[i]:
                neighbor3.append(list1[j])
            elif list1[j] == xlinked_list1[i]:
                if list2[j] != xlinked_list2[i]:
                    neighbor3.append(list2[j])
        for h in neighbor3:
            for l in range(0, len(list1)):
                if list1[l] == h:
                    if list2[l] != xlinked_list1[i]:
                        dihedrals.append([])
                        dihedrals[dc].append(xlinked_list2[i])
                        dihedrals[dc].append(xlinked_list1[i])
                        dihedrals[dc].append(h)
                        dihedrals[dc].append(list2[l])
                        dc += 1
                        neighbor4.append(list2[l])
                elif list2[l] == h:
                    if list1[l] != xlinked_list1[i]:
                        dihedrals.append([])
                        dihedrals[dc].append(xlinked_list2[i])
                        dihedrals[dc].append(xlinked_list1[i])
                        dihedrals[dc].append(h)
                        dihedrals[dc].append(list1[l])
                        dc += 1
                        neighbor4.append(list1[l])

        count = 0
        for m in neighbor4:
            pairs1.append(xlinked_list2[i])
            pairs2.append(neighbor4[count])
            #pairs.append('{:6d}{:7d}{:7d}'.format(list2[i], neighbor4[count], 1))
            count += 1
    return pairs1, pairs2, dihedrals


end1 = time.time()
print end1 - start

def uniq_pair(pairs1, pairs2, pairs_list_prev):
    for i in range(0, len(pairs1)):
        a = [pairs1[i], pairs2[i]]  # the order it is read from the pairs list
        b = [pairs2[i], pairs1[i]]
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

pairs1, pairs2, dihedrals = gen_pairs_list(list1, list2, xlinked_list1, xlinked_list2, nb)

unique_pairs1 = uniq_pair(pairs1, pairs2, pairs_list_prev)

#unique_dihedrals = uniq_dihedral(dihedrals, dihedrals_prev)
count = 0
# for i in unique_pairs:
#     count += 1
#
# print count

pairs1, pairs2, dihedrals = gen_pairs_list(list1, list2, xlinked_list3, xlinked_list4, nb)

unique_pairs2 = uniq_pair(pairs1, pairs2, pairs_list_prev)

unique_dihedrals = uniq_dihedral(dihedrals, dihedrals_prev)
count = 0
for i in dihedrals_prev:
    count += 1

print count

# for i in unique_pairs:
#     print i
#
# print end1 - start
# # end2 = time.time()
# # # print end2 - end1
# # # print end2 - start
# # # Now print them in order
# #
# # ordered_unique_pairs = []
# # for i in range(1, nr + 1):
# #     no_same_atom = []
# #     for k in range(0, len(unique_pairs)):
# #         if int(unique_pairs[k][0:6]) == i:
# #             no_same_atom.append(unique_pairs[k])
# #     for l in range(1, nr + 1):
# #         for j in range(0, len(no_same_atom)):
# #             if int(no_same_atom[j][6:13]) == l:
# #                 ordered_unique_pairs.append(no_same_atom[j])
# #
# # end3 = time.time()
# # print end3 - end2
# # print end3 - start
# #
# # # for i in range(0, len(ordered_unique_pairs)):
# # #     print ordered_unique_pairs[i]
# #
#print len(unique_pairs)
# #
# # TESTING
#
# pairs_index = 0  # find index where [ pairs ] section begins
# while a[pairs_index].count('[ pairs ]') == 0:
#     pairs_index += 1
#
# npair = 0  # number of lines in the 'pairs' section
# pairs_count = pairs_index + 2  # keep track of index of a
# while a[pairs_count] != '\n':
#     pairs_count += 1
#     npair += 1
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
