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

# pairs = []
# for i in range(bonds_index + 2, bond_count):  # look at each bond in [ bonds ] section
#     if_1 = 0
#     if_2 = 0
#     atom1 = int(a[i][0:6])  # Defining these reduces the time for this loop by about a third
#     atom2 = int(a[i][6:13])
#     for k in range(1, nr + 1):  # look for all atoms
#         neighbor3 = []  # 1-3 neighbors
#         neighbor4 = []  # 1-4 neighbors
#         if if_1 == 0:
#             if atom1 == k:
#                 if_1 += 1
#                 neighbor1 = atom2  # This the atom's immediate neighbor which it is bonded to
#                 for j in range(bonds_index + 2, bond_count):
#                     atom1 = a[j][0:6]
#                     atom2 = a[j][6:13]
#                     if atom1 == neighbor1:
#                         if atom2 != k:
#                             neighbor3.append(atom2)
#                     elif atom2 == neighbor1:
#                         if atom1 != k:
#                             neighbor3.append(atom1)
#                 for h in neighbor3:
#                     for l in range(bonds_index + 2, bond_count):
#                         atom1 = a[l][0:6]
#                         atom2 = a[l][6:13]
#                         if atom1 == h:
#                             if atom2 != neighbor1:
#                                 neighbor4.append(atom2)
#                         elif atom2 == h:
#                             if atom1 != neighbor1:
#                                 neighbor4.append(atom1)
#         if if_2 == 0:
#             if atom2 == k:
#                 if_2 += 1
#                 neighbor1 = atom1  # This the atom's immediate neighbor which it is bonded to
#                 for j in range(bonds_index + 2, bond_count):
#                     atom1 = a[j][0:6]
#                     atom2 = a[j][6:13]
#                     if atom1 == neighbor1:
#                         if atom2 != k:
#                             neighbor3.append(atom2)
#                     elif atom2 == neighbor1:
#                         if atom1 != k:
#                             neighbor3.append(atom1)
#                 for h in neighbor3:
#                     for l in range(bonds_index + 2, bond_count):
#                         atom1 = a[l][0:6]
#                         atom2 = a[l][6:13]
#                         if atom1 == h:
#                             if atom2 != neighbor1:
#                                 neighbor4.append(atom2)
#                         elif atom2 == h:
#                             if atom1 != neighbor1:
#                                 neighbor4.append(atom1)
#         count = 0
#         for m in neighbor4:
#             pairs.append('{:6d}{:7d}{:7d}'.format(k, neighbor4[count], 1))
#             count += 1
#         if if_1 == 1 and if_2 == 1:
#             break
start = time.time()

list1 = []  # list to hold atom numbers of 1st entry in each row of [ bonds ] section
list2 = []  # list to hold atom numbers of 2nd entry in each row of [ bonds ] section
for i in range(bonds_index + 2, bond_count):  # extract atom numbers in [ bonds ] section in order of appearance
    list1.append(int(a[i][0:6]))  # converts them to integers to save time later
    list2.append(int(a[i][6:13]))

max1 = max(list1)
max2 = max(list2)
max_overall = max(max1, max2)

# # Make adjacency matrix
# A = csr_matrix((65000, 65000), dtype=np.int8).toarray()
# A2 = A.dot(A)
# end = time.time()
# print end - start
# sys.exit()
# # for i in range(0, len(list1)):
# #     A[list1[i] - 1, list2[i] - 1] = 1
# #     A[list2[i] - 1, list1[i] - 1] = 1
# #
# # A2 = np.dot(A, A)
# # A3 = np.dot(A2, A)
#
# def gen_pairs(list1, list2):
#     max1 = max(list1)
#     max2 = max(list2)
#     max_overall = max(max1, max2)
#
#     # Make adjacency matrix
#     A = np.zeros((max_overall, max_overall))
#
#     for i in range(0, len(list1)):
#         A[list1[i] - 1, list2[i] - 1] = 1
#         A[list2[i] - 1, list1[i] - 1] = 1
#
#     A2 = np.dot(A, A)
#     A3 = np.dot(A2, A)
#     pairs = []
#     for i in range(0, max_overall):
#         for k in range(0, max_overall):
#             if A[i, k] == 0:
#                 if A2[i, k] == 0:
#                     if A3[i, k] != 0:
#                         pairs.append('{:6d}{:7d}{:7d}'.format(i + 1, k + 1, 1))
#     return pairs

# pairs = gen_pairs(list1, list2)

def gen_pairs_list(list1, list2):
    pairs = []
    for i in range(0, len(list1)):  # look at all bonds
        neighbor3 = []  # redefine lists for each loop since each loops looks at a different bond
        neighbor4 = []
        for j in range(0, len(list1)):  # compare to all other bonds
            if list1[j] == list2[i]:  # look for 1-3 pairs by comparing first entry of each row of [ pairs ] section
                neighbor3.append(list2[j])
            elif list2[j] == list2[i]:  # look for 1-3 pairs again
                if list1[j] != list1[i]:
                    neighbor3.append(list1[j])
        for h in neighbor3:
            for l in range(0, len(list1)):
                if list1[l] == h:
                    if list2[l] != list2[i]:
                        neighbor4.append(list2[l])
                elif list2[l] == h:
                    if list1[l] != list2[i]:
                        neighbor4.append(list1[l])

        count = 0
        for m in neighbor4:
            pairs.append('{:6d}{:7d}{:7d}'.format(list1[i], neighbor4[count], 1))
            count += 1

        neighbor3 = []
        neighbor4 = []
        for j in range(0, len(list1)):
            if list2[j] == list1[i]:
                neighbor3.append(list1[j])
            elif list1[j] == list1[i]:
                if list2[j] != list2[i]:
                    neighbor3.append(list2[j])
        for h in neighbor3:
            for l in range(0, len(list1)):
                if list1[l] == h:
                    if list2[l] != list1[i]:
                        neighbor4.append(list2[l])
                elif list2[l] == h:
                    if list1[l] != list1[i]:
                        neighbor4.append(list1[l])

        count = 0
        for m in neighbor4:
            pairs.append('{:6d}{:7d}{:7d}'.format(list2[i], neighbor4[count], 1))
            count += 1
    return pairs

pairs = gen_pairs_list(list1, list2)

# print pairs
# end1 = time.time()
#
# print end1 - start
#
# # eliminate duplicates
#
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

for i in unique_pairs:
    print i

# end2 = time.time()
# # print end2 - end1
# # print end2 - start
# # Now print them in order
#
# ordered_unique_pairs = []
# for i in range(1, nr + 1):
#     no_same_atom = []
#     for k in range(0, len(unique_pairs)):
#         if int(unique_pairs[k][0:6]) == i:
#             no_same_atom.append(unique_pairs[k])
#     for l in range(1, nr + 1):
#         for j in range(0, len(no_same_atom)):
#             if int(no_same_atom[j][6:13]) == l:
#                 ordered_unique_pairs.append(no_same_atom[j])
#
# end3 = time.time()
# print end3 - end2
# print end3 - start
#
# # for i in range(0, len(ordered_unique_pairs)):
# #     print ordered_unique_pairs[i]
#
print len(unique_pairs)
# print end1 - start
#
# TESTING

pairs_index = 0  # find index where [ pairs ] section begins
while a[pairs_index].count('[ pairs ]') == 0:
    pairs_index += 1

npair = 0  # number of lines in the 'pairs' section
pairs_count = pairs_index + 2  # keep track of index of a
while a[pairs_count] != '\n':
    pairs_count += 1
    npair += 1

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

if count == 349:
    print 'It works!'
