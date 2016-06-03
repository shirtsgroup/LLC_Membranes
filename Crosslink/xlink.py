# Script to crosslink LLC monomers based on distance between carbons at the end of a simulation

#!/usr/bin/python

import os
import argparse
import numpy as np
import math
import matplotlib.pyplot as plt

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')

# Flags
parser.add_argument('-i', '--input', default='wiggle.gro', help = 'Name of input file')
parser.add_argument('-t', '--type', default='LLC', help='Type of monomer')
parser.add_argument('-l', '--layers', default=20, help='Number of layers')
parser.add_argument('-p', '--pores', default=4, help='Number of pores')
parser.add_argument('-n', '--no_monomers', default=6, help='Number of monomers per layer')
parser.add_argument('-c', '--cutoff', default=0.313, help='Anything below this value will be cross linked')

args = parser.parse_args()

# Read in coordinate file from simulation

f = open('wiggle.gro', 'r')
a = []
lines = 0
for line in f:
    a.append(line)
    lines += 1  # count number of lines in the file

lines_of_text = 0
for i in range(0, len(a)):
    if a[i].count('1LLC') == 0:
        lines_of_text += 1
    if a[i].count('1LLC') == 1:
        break

if args.type == 'LLC':
    atoms = 138
    atoms_list = ['O5', 'O6', 'O7', '08', 'O9', 'C32', 'C18', 'C46', 'C33', 'C19', 'C47', 'C34', 'C20', 'C47', 'H49',
                  'H24', 'H74', 'H50', 'H25', 'H75', 'H51', 'H26', 'H76']
    topology = 'HII_mon.itp'

# Use Assembly_itp.py to create a topology for the entire assembly and then write it to a file which will be edited to
# incorporate cross-links

import subprocess
with open("crosslinked.itp", "w+") as output:
    subprocess.call(["python", "./../Structure-Files/Assembly_itp.py"], stdout=output);

tot_atoms = atoms*args.layers*args.no_monomers*args.pores
tot_monomers = args.layers*args.no_monomers*args.pores

# List of coordinates for carbon 1 (i.e. the carbon on the end of the chain)
C1x = []  # list to hold x positions of atoms of interest
C1y = []  # list to hold y positions of atoms of interest
C1z = []  # list to hold z positions of atoms of interest. Z axis runs parallel to pore

# List of coordinates for carbon 2 (i.e. the second carbon from the end of the chains
C2x = []  # list to hold x positions of atoms of interest
C2y = []  # list to hold y positions of atoms of interest
C2z = []  # list to hold z positions of atoms of interest. Z axis runs parallel to pore

for line in range(0, lines):  # 15 is where the field containing atom identities ends
    if a[line][0:15].count('C20') == 1 or a[line][0:15].count('C34') == 1 or a[line][0:15].count('C48') == 1:
        C1x.append(float(a[line][20:28]))
        C1y.append(float(a[line][28:36]))
        C1z.append(float(a[line][36:44]))
    elif a[line][0:15].count('C19') == 1 or a[line][0:15].count('C33') == 1 or a[line][0:15].count('C47') == 1:
        C2x.append(float(a[line][20:28]))
        C2y.append(float(a[line][28:36]))
        C2z.append(float(a[line][36:44]))

# Find distance between carbon 1 and carbon 2 for all pairs

dist = np.zeros((len(C1x), len(C2x)))

for i in range(0, len(C1x)):
    for k in range(0, len(C2x)):
        if i != k:  # make sure that C1 and C2 that are a part of the same monomer do not factor into this calculation
            dist[k, i] = math.sqrt((C1x[i] - C2x[k])**2 + (C1y[i] - C2y[k])**2 + (C1z[i] - C2z[k])**2)
        else:
            dist[k, i] = 4  # artificially high number to keep it from interfering

# Find the distance of the closest carbon and its index
min_dist = np.zeros((len(C1x), 1))
min_index = np.zeros((len(C1x), 1))  # index of minimum value of distances for each monomer-monomer measurement
for i in range(0, len(C1x)):
    min_dist[i, 0] = min(dist[:, i])
    min_index[i, 0] = np.argmin(dist[:, i])  # index corresponds to the monomer with which the minimum C1-C2 distance is achieved

# Now see which of these distances meet the cutoff criteria

change_dist = []  # distances associated with atoms which need to be changed
change_index1 = []  # index of C1 from min_dist which needs to be changed
change_index2 = []  # index of C2 from dist which needs to be changed (index already contained in min_index but will be
                    # truncated here

count = 0
for i in range(0, len(C1x)):
    if min_dist[i, 0] <= args.cutoff:
        change_dist.append(min_dist[i, 0])
        change_index1.append(i)
        change_index2.append(int(min_index[i, 0]))
        count += 1

# Now that everything has an index that needs to be changed, we must interpret those indices

# C19 is the 26th atom, C20 is 27th, C33 is 41st, C34 is 42nd, C47 is 56th, C48 is 57th, in each monomer

# We also need to know which index refers to which monomer and tail -- applies equally for C1 and C2
# Every third index starts a new monomer. Each index in between is a tail (inclusive)

tail1 = np.arange(0, len(C1x) + 1, 3)  # indices of carbons
tail2 = np.arange(1, len(C1x) + 1, 3)
tail3 = np.arange(2, len(C1x) + 1, 3)

C1_no = []
C2_no = []
for i in range(0, count):
    if change_index1[i] in tail1:  # This is C1 therefore if this is true, then the atom is C20 (atom no 27)
        C1_no.append((change_index1[i]/3)*(atoms - 1) + 27)
    if change_index1[i] in tail2:  # This is C1 therefore if this is true, then the atom is C34 (atom no 42)
        C1_no.append((change_index1[i]/3)*(atoms - 1) + 42)  # division automatically round down
    if change_index1[i] in tail3:  # This is C1 therefore if this is true, then the atom is C48 (atom no 57)
        C1_no.append((change_index1[i]/3)*(atoms - 1) + 57)  # division automatically round down
    if change_index2[i] in tail1:  # This is C2 therefore if this is true, then the atom is C19 (atom no 26)
        C2_no.append((change_index2[i]/3)*(atoms - 1) + 26)
    if change_index2[i] in tail2:  # This is C2 therefore if this is true, then the atom is C33 (atom no 41)
        C2_no.append((change_index2[i]/3)*(atoms - 1) + 41)  # division automatically round down
    if change_index2[i] in tail3:  # This is C2 therefore if this is true, then the atom is C47 (atom no 56)
        C2_no.append((change_index2[i]/3)*(atoms - 1) + 56)  # division automatically round down

print C1_no, C2_no
print change_dist
print max(dist[0, :])