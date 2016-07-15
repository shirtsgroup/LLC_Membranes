# Script to crosslink LLC monomers based on distance between carbons at the end of a simulation

#!/usr/bin/python

import os
import argparse
import numpy as np
import math
import time
import genpairs
import copy

start = time.time()  # For informational purposes
location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')  # allow input from user

# Flags
parser.add_argument('-i', '--input', default='wiggle_init.gro', help = 'Name of input file')
parser.add_argument('-t', '--type', default='LLC', help='Type of monomer')
parser.add_argument('-l', '--layers', default=20, help='Number of layers')
parser.add_argument('-p', '--pores', default=4, help='Number of pores')
parser.add_argument('-n', '--no_monomers', default=6, help='Number of monomers per layer')
parser.add_argument('-c', '--cutoff', default=5, help='Cutoff distance for cross-linking. Bottom x % of the distribution'
                                                      'of distances will be cross-linked ')
parser.add_argument('-e', '--term_prob', default=5, help='Termination probability (%)')

args = parser.parse_args()

# Read in coordinate file from simulation

f = open(args.input, 'r')  # open the file which has just finished being simulated and copy in all of the information
a = []
lines = 0
for line in f:
    a.append(line)
    lines += 1  # count number of lines in the file

lines_of_text = 0  # find the number of lines of relevant text
for i in range(0, len(a)):
    if a[i].count('1LLC') == 0:
        lines_of_text += 1
    if a[i].count('1LLC') == 1:
        break

if args.type == 'LLC':
    atoms = 143  # not including sodium but including dummy atoms
    atoms_list = ['O5', 'O6', 'O7', '08', 'O9', 'C32', 'C18', 'C46', 'C33', 'C19', 'C47', 'C34', 'C20', 'C47', 'H49',
                  'H24', 'H74', 'H50', 'H25', 'H75', 'H51', 'H26', 'H76']
    topology = 'HII_mon.itp'
    xlink_atoms = 6  # number of atoms involved in cross-linking
    no_dummies = 6


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

no_carbons = len(C1x)  # number of carbon atoms that may be involved in cross-linking

# Time to handle periodic boundary conditions

x1 = float(a[len(a) - 1][0:10])  # x box vector
y1 = float(a[len(a) - 1][10:20])  # component one of y box vector
y2 = float(a[len(a) - 1][50:60])  # component two of y box vector

# Visualize pbc's like this for the following calculations:
#  ___________________________
#  \        \        \        \
#   \   I    \   II   \  III   \
#    \________\________\________\
#     \        \        \        \
#      \  IV    \  Main  \   V    \
#       \________\________\________\
#        \        \        \        \
#         \   VI   \  VII   \  VIII  \
#          \________\________\________\
#
# We will use x1, y1 and y2 to make copies of the coordinates in the 'Main' box into each of the labeled boxes

# Define lists to hold translated coordinates
# Lists will be of the format [C1x C1y C1z C2x C2y C2z]
I = []
II = []
III = []
IV = []
V = []
VI = []
VII = []
VIII = []

for i in range(0, xlink_atoms):
    I.append([]), II.append([]), III.append([]), IV.append([]), V.append([]), VI.append([]), VII.append([]), VIII.append([])

for i in range(0, len(C1x)):
    I[0].append(C1x[i] - x1 + y2), I[1].append(C1y[i] + y1), I[2].append(C1z[i]),\
    I[3].append(C2x[i] - x1 + y2), I[4].append(C2y[i] + y1), I[5].append(C2z[i])
    II[0].append(C1x[i] + y2), II[1].append(C1y[i] + y1), II[2].append(C1z[i]),\
    II[3].append(C2x[i] + y2), II[4].append(C2y[i] + y1), II[5].append(C2z[i])
    III[0].append(C1x[i] + x1 + y2), III[1].append(C1y[i] + y1), III[2].append(C1z[i]), \
    III[3].append(C2x[i] + x1 + y2), III[4].append(C2y[i] + y1), III[5].append(C2z[i])
    IV[0].append(C1x[i] - x1), IV[1].append(C1y[i]), IV[2].append(C1z[i]),\
    IV[3].append(C2x[i] - x1), IV[4].append(C2y[i]), IV[5].append(C2z[i])
    V[0].append(C1x[i] + x1), V[1].append(C1y[i]), V[2].append(C1z[i]),\
    V[3].append(C2x[i] + x1), V[4].append(C2y[i]), V[5].append(C2z[i])
    VI[0].append(C1x[i] - x1 - y2), VI[1].append(C1y[i] - y1), VI[2].append(C1z[i]),\
    VI[3].append(C2x[i] - x1 - y2), VI[4].append(C2y[i] - y1), VI[5].append(C2z[i])
    VII[0].append(C1x[i] - y2), VII[1].append(C1y[i] - y1), VII[2].append(C1z[i]), \
    VII[3].append(C2x[i] - y2), VII[4].append(C2y[i] - y1), VII[5].append(C2z[i])
    VIII[0].append(C1x[i] + x1 - y2), VIII[1].append(C1y[i] - y1), VIII[2].append(C1z[i]), \
    VIII[3].append(C2x[i] + x1 - y2), VIII[4].append(C2y[i] - y1), VIII[5].append(C2z[i])

# Now all of the translated coordinates need to be added to the main list of x, y and z coordinates for each carbon

for i in range(0, len(I[0])):
    C1x.append(I[0][i]), C1y.append(I[1][i]), C1z.append(I[2][i]), C2x.append(I[3][i]), C2y.append(I[4][i])
    C2z.append(I[5][i])

for i in range(0, len(II[0])):
    C1x.append(II[0][i]), C1y.append(II[1][i]), C1z.append(II[2][i]), C2x.append(II[3][i]), C2y.append(II[4][i])
    C2z.append(II[5][i])

for i in range(0, len(III[0])):
    C1x.append(III[0][i]), C1y.append(III[1][i]), C1z.append(III[2][i]), C2x.append(III[3][i]), C2y.append(III[4][i])
    C2z.append(III[5][i])

for i in range(0, len(IV[0])):
    C1x.append(IV[0][i]), C1y.append(IV[1][i]), C1z.append(IV[2][i]), C2x.append(IV[3][i]), C2y.append(IV[4][i])
    C2z.append(IV[5][i])

for i in range(0, len(V[0])):
    C1x.append(V[0][i]), C1y.append(V[1][i]), C1z.append(V[2][i]), C2x.append(V[3][i]), C2y.append(V[4][i])
    C2z.append(V[5][i])

for i in range(0, len(VI[0])):
    C1x.append(VI[0][i]), C1y.append(VI[1][i]), C1z.append(VI[2][i]), C2x.append(VI[3][i]), C2y.append(VI[4][i])
    C2z.append(VI[5][i])

for i in range(0, len(VII[0])):
    C1x.append(VII[0][i]), C1y.append(VII[1][i]), C1z.append(VII[2][i]), C2x.append(VII[3][i]), C2y.append(VII[4][i])
    C2z.append(VII[5][i])

for i in range(0, len(VIII[0])):
    C1x.append(VIII[0][i]), C1y.append(VIII[1][i]), C1z.append(VIII[2][i]), C2x.append(VIII[3][i]), C2y.append(VIII[4][i])
    C2z.append(VIII[5][i])
stop1 = time.time()
print 'PBCs set up: %s seconds' %(stop1 - start)

# Find distance between carbon 1 and carbon 2 for all pairs

dist = np.zeros((len(C1x), len(C2x)))

# Create a matrix of 1's and 0's. Entries that are 1's should not be counted in the distance calculations because they
# represent distances that are between carbons and the same monomer

a = np.ones((1, 12960))[0]
b = np.ones((1, 11520))[0]
c = np.ones((1, 10080))[0]
d = np.ones((1, 8640))[0]
e = np.ones((1, 7200))[0]
f = np.ones((1, 5760))[0]
g = np.ones((1, 4320))[0]
h = np.ones((1, 2880))[0]
i = np.ones((1, 1440))[0]

exclude = np.diag(a, 0) + np.diag(b, -1440) + np.diag(b, 1440) + np.diag(c, 2880) + np.diag(c, -2880) + np.diag(d, 4320) + \
    np.diag(d, -4320) + np.diag(e, 5760) + np.diag(e, -5760) + np.diag(f, 7200) + np.diag(f, -7200) + \
    np.diag(g, 8640) + np.diag(g, -8640) + np.diag(h, 10080) + np.diag(h, -10080) + np.diag(i, 11520) + np.diag(i, -11520)

# for i in range(0, len(C1x)):
#     for k in range(0, len(C2x)):
#         if exclude[k, i] == 1:  # make sure that C1 and C2 that are a part of the same monomer do not factor into this calculation
#             dist[k, i] = 1000  # artificially high number to keep it from interfering
#         else:
#             dist[k, i] = math.sqrt((C1x[i] - C2x[k])**2 + (C1y[i] - C2y[k])**2 + (C1z[i] - C2z[k])**2)

dist = genpairs.calc_dist(C1x, C2x, C1y, C2y, C1z, C2z, exclude, dist)

stop2 = time.time()

print 'Distances Calculated: %s seconds' %(stop2 - stop1)

# Find the distance of the closest carbon and its index
min_dist = np.zeros((len(C1x), 1))
min_index = np.zeros((len(C1x), 1))  # index of minimum value of distances for each monomer-monomer measurement
for i in range(0, len(C1x)):
    min_dist[i, 0] = min(dist[:, i])
    min_index[i, 0] = np.argmin(dist[:, i]) % (len(C1x)/9)  # index corresponds to the monomer with which
                                                                        # the minimum C1-C2 distance is achieved

# Now see which of these distances meet the cutoff criteria

change_index1 = []  # index of C1 from min_dist which needs to be changed
change_index2 = []  # index of C2 from dist which needs to be changed (index already contained in min_index but will be
                    # truncated here)

# find the cutoff distance for cross-linking

min_list = min_dist.tolist()
for i in range(0, int((float(args.cutoff)/100)*len(C1x))):  # looks at a percentage of the total values based on user input of cutoff
    m = min(min_list)  # finds minimum of min_list
    min_list.remove(m)  # removes that value from min_list

cutoff = min(min_list)  # The minimum value left after the modification of min_list is the cutoff value

for i in range(0, len(C1x)):
    if min_dist[i, 0] < cutoff:  # find distances which meet the cutoff criteria
        change_index1.append(i % (len(C1x)/9))  # holds the index of the atom associated with the met criteria for primary carbons (C20, C34, C48)
        change_index2.append(int(min_index[i, 0]))  # Same as above but for the secondary carbons

# Now remove duplicates while preserving order (if there are any)

# change_index1 = uniq(change_index1)
# change_index2 = uniq(change_index2)

# Now that everything has an index that needs to be changed, we must interpret those indices

# C19 is the 26th atom, C20 is 27th, C33 is 41st, C34 is 42nd, C47 is 56th, C48 is 57th, in each monomer

# We also need to know which index refers to which monomer and tail -- applies equally for C1 and C2
# Every third index starts a new monomer. Each index in between is a tail (inclusive)

tail1 = np.arange(0, len(C1x) + 1, 3)  # indices of carbons
tail2 = np.arange(1, len(C1x) + 1, 3)
tail3 = np.arange(2, len(C1x) + 1, 3)

C1_no = []
C2_no = []
for i in range(0, len(change_index1)):
    if change_index1[i] in tail1:  # This is C1 therefore if this is true, then the atom is C20 (atom no 27)
        C1_no.append((change_index1[i]/3)*atoms + 27)
    if change_index1[i] in tail2:  # This is C1 therefore if this is true, then the atom is C34 (atom no 42)
        C1_no.append((change_index1[i]/3)*atoms + 42)  # division automatically round down
    if change_index1[i] in tail3:  # This is C1 therefore if this is true, then the atom is C48 (atom no 57)
        C1_no.append((change_index1[i]/3)*atoms + 57)  # division automatically round down
    if change_index2[i] in tail1:  # This is C2 therefore if this is true, then the atom is C19 (atom no 26)
        C2_no.append((change_index2[i]/3)*atoms + 26)
    if change_index2[i] in tail2:  # This is C2 therefore if this is true, then the atom is C33 (atom no 41)
        C2_no.append((change_index2[i]/3)*atoms + 41)  # division automatically round down
    if change_index2[i] in tail3:  # This is C2 therefore if this is true, then the atom is C47 (atom no 56)
        C2_no.append((change_index2[i]/3)*atoms + 56)  # division automatically round down

# reduce lists down to unique atom pairs


def uniq(input1, input2):
    output1 = []
    output2 = []
    for i in range(0, len(input1)):
        if input1[i] not in output1 and input2[i] not in output2:
            output1.append(input1[i])
            output2.append(input2[i])
    return output1, output2

c1_uniq, c2_uniq = uniq(C1_no, C2_no)

# In some cases, it is possible that a c1 and a c2 on the same chain may be involved in independent cross-link reactions
# We need to get rid of those to avoid unnecessary confusion. Since c2 is the reactive atom (it is a radical at some
# point), it will trump c1. So if a c2 and c1 are on the same chain, then c1 will be removed from consideration.

# Create a list of the indices of all the c1 atoms adjacent to reacting c2's
adjacent_c1 = []
for i in range(0, len(c2_uniq)):
    adjacent_c1.append(c2_uniq[i] + 1)

# Create the final lists of c1 and c2 that will be cross-linked
c1 = []
c2 = []
for i in range(0, len(adjacent_c1)):
    if c1_uniq[i] not in adjacent_c1:  # add values to the lists that have c1's that are not adjacent to reacting c2's
        c1.append(c1_uniq[i])
        c2.append(c2_uniq[i])

# Some of the reactions will terminate instead of cross linking. This will result in the dummy hydrogen on C2 to become
# active.

tp = args.term_prob  # termination probability
if type(tp) == int:
    no_decimals = 0  # there aren't decimals in a integer
else:
    no_decimals = len(str(tp)) - len(str(int(tp))) - 1  # subtract the length of the integer value from the length
                                                        # of the float value and subtract one for the decimal point to
                                                        # get the number of decimal places in the number

# We must make an appropriate length array of 0's and 1's in order to create a probability distribution
# If tp is 5 %, then it should be an array of length 20, with 19 0's and a 1
# If tp is 5.46 % then it should be an array of length 10000 with 546 1's and the rest 0

array_length = 100 * 10 ** no_decimals  # multiply 100 % by 10 raised to the number of decimal points
no_ones = int(tp * 10 ** no_decimals)  # multiply the termination probability by the same amount

term_prob_array = []  # make an array with a a percentage of 1's equal to the termination probability
for i in range(0, array_length - no_ones):
    term_prob_array.append(0)

for i in range(array_length - no_ones, array_length):
    term_prob_array.append(1)

# Termination occurs at C2
# Move termination carbons to a new list
term_c1 = []
term_c2 = []
term_indices = []
no_xlinks = copy.deepcopy(len(c2))  # The length of c2 may change in the next loops so I am preserving it here

# Add the termination carbons to new list
for i in range(0, no_xlinks):
    term = np.random.choice(term_prob_array)
    if term == 1:
        term_c1.append(c1[i])
        term_c2.append(c2[i])
        term_indices.append(i)

# Remove terminated carbons from c1 and c2
for i in range(0, len(term_indices)):
    len_diff = no_xlinks - len(c1)  # Account for varying size of c1 and c2 lists
    del c1[term_indices[i] - len_diff]
    del c2[term_indices[i] - len_diff]

# Use Assembly_itp.py to create a topology for the entire assembly and then write it to a file which will be edited to
# incorporate cross-links

import subprocess
with open("crosslinked.itp", "w+") as output:
    subprocess.call(["python", "./../Structure-Files/Assembly_itp.py", "-x", "on"], stdout=output);

# open and read that new file

f = open('crosslinked.itp', 'r')
print 'Crosslinked.itp file written'
b = []
for line in f:
    b.append(line)

# Make lists with numbers of H atoms whose type needs to change as a consequence of the cross linking reaction
# Type should change from ha to hc

# first, redefine tail 1, 2 and 3 using C2 atom numbers as a reference
# find the residue numbers of the first C33, C47, C19 (C2 carbons)

index = 0
while str.strip(b[index][22:29]) != 'C19':
    index += 1  # increments the while loop

C2_1 = int(b[index][0:5])

index = 0
while str.strip(b[index][22:29]) != 'C33':
    index += 1  # increments the while loop

C2_2 = int(b[index][0:5])

index = 0
while str.strip(b[index][22:29]) != 'C47':
    index += 1  # increments the while loop

C2_3 = int(b[index][0:5])

tail1 = []
tail2 = []
tail3 = []

for i in range(0, tot_monomers):
    tail1.append(atoms*i + C2_1)
    tail1.append(atoms*i + C2_1 + 1)
    tail2.append(atoms*i + C2_2)
    tail2.append(atoms*i + C2_2 + 1)
    tail3.append(atoms*i + C2_3)
    tail3.append(atoms*i + C2_3 + 1)

H = []  # list of residue numbers of H that need to be changed from ha to hc
H_new1 = []  # list of residue number of H's that will be convert from dummy atoms to real atoms
for i in range(0, len(c2)):
    if c2[i] in tail1:  # i.e. C19
        H.append(c2[i] + 59)
        H.append(c2[i] + 60)
        H.append(c2[i] + 61)
        if c1[i] in tail1:
            H.append(c1[i] + 59)
            H.append(c1[i] + 60)
        if c1[i] in tail2:
            H.append(c1[i] + 69)
            H.append(c1[i] + 70)
        if c1[i] in tail3:
            H.append(c1[i] + 79)
            H.append(c1[i] + 80)
        H_new1.append(c2[i] + 113)
    if c2[i] in tail2:  # i.e. C33
        H.append(c2[i] + 69)
        H.append(c2[i] + 70)
        H.append(c2[i] + 71)
        if c1[i] in tail1:
            H.append(c1[i] + 59)
            H.append(c1[i] + 60)
        if c1[i] in tail2:
            H.append(c1[i] + 69)
            H.append(c1[i] + 70)
        if c1[i] in tail3:
            H.append(c1[i] + 79)
            H.append(c1[i] + 80)
        H_new1.append(c2[i] + 100)
    if c2[i] in tail3:  # i.e. C47
        H.append(c2[i] + 79)
        H.append(c2[i] + 80)
        H.append(c2[i] + 81)
        if c1[i] in tail1:
            H.append(c1[i] + 59)
            H.append(c1[i] + 60)
        if c1[i] in tail2:
            H.append(c1[i] + 69)
            H.append(c1[i] + 70)
        if c1[i] in tail3:
            H.append(c1[i] + 79)
            H.append(c1[i] + 80)
        H_new1.append(c2[i] + 87)

# need to do the same for the termination carbons except we need to include the dummy hydrogen and don't have to worry
# about the term_c1 list since nothing is happening on that end anymore

# the c2 where termination is happening is no longer the same c2 which is referenced in term_c2
term = []  # define those termination atoms as 'term' which are the carbons adjacent to term_c1
for i in range(0, len(term_c1)):
    term.append(term_c1[i] - 1)

H_new2 = []  # keep these separate for indexing purposes
for i in range(0, len(term)):
    if term[i] in tail1:  # i.e. C19
        H.append(term[i] + 59)
        H.append(term[i] + 60)
        H.append(term[i] + 61)
        H_new2.append(term[i] + 112)
        H_new2.append(term[i] + 113)
    if term[i] in tail2:  # i.e. C33
        H.append(term[i] + 69)
        H.append(term[i] + 70)
        H.append(term[i] + 71)
        H_new2.append(term[i] + 99)
        H_new2.append(term[i] + 100)
    if term[i] in tail3:  # i.e. C47
        H.append(term[i] + 79)
        H.append(term[i] + 80)
        H.append(term[i] + 81)
        H_new2.append(term[i] + 86)
        H_new2.append(term[i] + 87)

# Before going further, we need to label all of the dummy atoms that are still dummies so that they are not written into
# the bonds, angles, dihedrals or pairs section

# generate a list of indices of the dummy atoms
dummy_list = []
for i in range(0, tot_monomers):
    for k in range(0, no_dummies):
        dummy_list.append(i*atoms + (atoms - k))  # dummies are listed at the end of the atomtypes

leftovers = []  # poor dummies that didn't get turned into real atoms
for i in range(0, len(dummy_list)):
    if dummy_list[i] not in H_new1 and dummy_list[i] not in H_new2:  # sees if dummies were turned to real atoms
        leftovers.append(dummy_list[i])  # if not, they become part of the leftovers

# find the indices of all fields that need to be modified

# [ atoms ]

atoms_index = 0  # find index where [ atoms ] section begins
while b[atoms_index].count('[ atoms ]') == 0:
    atoms_index += 1

atoms_count = atoms_index + 2
nr = 0  # number of lines in 'atoms' section
while b[atoms_count] != '\n':
    atoms_count += 1  # increments the while loop
    nr += 1  # counts number of atoms

# [ bonds ]

bonds_index = 0  # find index where [ bonds ] section begins\
while b[bonds_index].count('[ bonds ]') == 0:
    bonds_index += 1

nb = 0  # number of lines in the 'bonds' section
bond_list = []
bond_count = bonds_index + 2
while b[bond_count] != '\n':
    a_b = int(b[bond_count][0:6])
    b_b = int(b[bond_count][6:13])
    if a_b not in leftovers and b_b not in leftovers:
        bond_list.append([a_b, b_b])
        nb += 1  # counting number of lines in 'bonds' section
    bond_count += 1  # increments while loop

# properly format bonds list
bonds = []
for i in range(0, len(bond_list)):
    bonds.append('{:6d}{:7d}{:4d}'.format(bond_list[i][0], bond_list[i][1], 1))

# eliminate the old bonds list

del b[bonds_index + 2: bond_count]

# Insert new bonds list

count = bonds_index + 2
for i in range(0, len(bonds)):
    b.insert(count, bonds[i] + '\n')
    count += 1

# Replace atom types of bonding carbons with sp3 hybridized carbons

# To be safe, we also need to change the atom type of the non-bonding c1. If it is initiated, then it goes from c2 to c3
other_c1 = []
radical_c2 = []
for i in range(0, len(c2)):
    other_c1.append(c2[i] + 1)
    radical_c2.append(c1[i] - 1)

# The other c2 keeps its hybridization

# Now make all of the changes
# Also, any carbons that are termination sites and/or become sp3 hybridized will be marked with a 'T' since it should
# not be considered in further cross-link iterations
count = 0
for i in range(atoms_index + 2, atoms_count):
    res_num = int(b[i][0:5])
    if res_num in c1:  # change bonding carbon 1 from c2 to c3 since it becomes sp3 hybridized
        b[i] = b[i][0:5] + b[i][5:10].replace('c2', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + 'T' + '\n'
    if res_num in other_c1:  # The c1 on the same tail as the bonding c2 is now sp3 hybridized
        b[i] = b[i][0:5] + b[i][5:10].replace('c2', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + 'T' + '\n'
    if res_num in c2:  # the bonding c2 becomes sp3 hybridized
        b[i] = b[i][0:5] + b[i][5:10].replace('ce', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + 'T' + '\n'
    if res_num in term:  # the terminating c2 (where the termination actually happens) becomes sp3 hybridized
        b[i] = b[i][0:5] + b[i][5:10].replace('ce', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + 'T' + '\n'
    # if res_num in H:  # all of the H's that attach to the new sp3 hybridized carbons
    #     b[i] = b[i][0:5] + b[i][5:10].replace('ha', 'hc') + b[i][10:len(b[atoms_index + 2])]
    if res_num in H_new1:  # Dummy atoms become real during the bonding process and have mass
        b[i] = b[i][0:53] + b[i][53:61].replace('0.00000', '1.00800') + b[i][61:len(b[atoms_index + 2])]
    if res_num in H_new2:  # Dummy atoms which become real during the termination process and have mass
        b[i] = b[i][0:53] + b[i][53:61].replace('0.00000', '1.00800') + b[i][61:len(b[atoms_index + 2])]
    if res_num in radical_c2:  # mark the c2 atoms that are now radicals. They remain sp2 hybridized but are reactive
        b[i] = b[i][0:5] + b[i][5:10].replace('ce', 'c2') + b[i][10:len(b[atoms_index + 2]) - 1] + '*' + '\n'
    if res_num in term_c1:  # c1 in a terminating chain will be sp3
        b[i] = b[i][0:5] + b[i][5:10].replace('c2', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + 'T' + '\n'

# Add bonds between cross-linking atoms

for i in range(0, len(c1)):
    b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(c1[i], c2[i], 1) + "\n")
    nb += 1

# Add new H bonds formed during propagation (dummy H's becoming real)
for i in range(0, len(H_new1)):
    b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(c1[i], H_new1[i], 1) + "\n")
    nb += 1

# new H bonds formed during termination (2 dummy H's becoming real)
k = 0
for i in range(0, len(term)):
    b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(term[i], H_new2[k], 1) + "\n")
    k += 1
    nb += 1
    b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(term[i] + 1, H_new2[k], 1) + "\n")
    nb += 1
    k += 1

# Now we need to generate dihedrals, angles and pairs lists. We have all the atoms we need. Let's make one list that
# holds all of these values:

atoms_of_interest = radical_c2 + c1 + c2 + other_c1 + term + term_c1 + H + H_new1 + H_new2

# See genpairs.pyx to see what is going on in the following lines:

# Generate a list of all bonds. 'list' contains all bonds. 'Condensed_list' contains bonds only relevant to 'atoms_of_interest'
list, condensed_list = genpairs.read_lists(b, atoms_of_interest, bonds_index + 2, nb)

# Generate a pairs list, a dihedrals list, and an angles list that includes those formed by the newly created bonds

pairs_list, dihedrals_list, angles_list = genpairs.gen_pairs_list(list, condensed_list)

# Find where the [ pairs ] section is located

pairs_index = 0  # find index where [ pairs ] section begins
while b[pairs_index].count('[ pairs ]') == 0:
    pairs_index += 1

npair = 0  # number of lines in the 'pairs' section
pairs_count = pairs_index + 2  # keep track of index of a
pairs_list_prev = []
while b[pairs_count] != '\n':
    a_p = b[pairs_count][0:6]
    b_p = b[pairs_count][6:13]
    if a_p not in leftovers and b_p not in leftovers:
        pairs_list_prev.append([a_p, b_p])  # record pairs that already exist
        npair += 1
    pairs_count += 1

# Eliminate duplicates in the entire pairs list
unique_pairs = genpairs.uniq_pair(pairs_list, pairs_list_prev)

# Create a properly formatted pairs list
pairs = []
for i in range(0, len(unique_pairs)):
    pairs.append('{:6s}{:7s}{:7d}'.format(unique_pairs[i][0], unique_pairs[i][1], 1))

# eliminate the old pairs list

del b[pairs_index + 2: pairs_count]

# Insert new pairs list

count = pairs_index + 2
for i in range(0, len(pairs)):
    b.insert(count, pairs[i] + '\n')
    count += 1

# [ angles ]

angles_index = 0  # find index where [ angles ] section begins
while b[angles_index].count('[ angles ]') == 0:
    angles_index += 1

na = 0  # number of lines in the 'angles' section
angle_count = angles_index + 2  # keep track of index of a
angles_prev = []
while b[angle_count] != '\n':
    a_a = int(b[angle_count][0:6])
    b_a = int(b[angle_count][6:13])
    c_a = int(b[angle_count][13:20])
    if a_a not in leftovers and b_a not in leftovers and c_a not in leftovers:
        angles_prev.append([a_a, b_a, c_a])
        na += 1
    angle_count += 1

# Eliminate duplicates in the entire angles list
unique_angles = genpairs.uniq_angles(angles_list, angles_prev)

# Create a properly formatted angles list
angles = []
for i in range(0, len(unique_angles)):
    angles.append('{:6d}{:7d}{:7d}{:7d}'.format(unique_angles[i][0], unique_angles[i][1],
                                                        unique_angles[i][2], 1))

# eliminate the old angles list

del b[angles_index + 2: angle_count]

# Insert new angles list

count = angles_index + 2
for i in range(0, len(angles)):
    b.insert(count, angles[i] + "\n")
    count += 1

# [ dihedrals ] ; propers

dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
while b[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
    dihedrals_p_index += 1

ndp = 0  # number of lines in the 'dihedrals ; proper' section
dihedrals_p_count = dihedrals_p_index + 2  # keep track of index of a
dihedrals_prev = []
while b[dihedrals_p_count] != '\n':
    a_d = int(b[dihedrals_p_count][0:6])
    b_d = int(b[dihedrals_p_count][6:13])
    c_d = int(b[dihedrals_p_count][13:20])
    d_d = int(b[dihedrals_p_count][20:27])
    if a_d not in leftovers and b_d not in leftovers and c_d not in leftovers and d_d not in leftovers:
        dihedrals_prev.append([a_d, b_d, c_d, d_d])
        ndp += 1
    dihedrals_p_count += 1

# Eliminate duplicates in entire dihedrals list
unique_dihedrals = genpairs.uniq_dihedral(dihedrals_list, dihedrals_prev)

# properly format dihedrals list
dihedrals = []
for i in range(0, len(unique_dihedrals)):
    dihedrals.append('{:6d}{:7d}{:7d}{:7d}{:4d}'.format(unique_dihedrals[i][0], unique_dihedrals[i][1],
                                                        unique_dihedrals[i][2], unique_dihedrals[i][3], 3))

# eliminate the old dihedrals list

del b[dihedrals_p_index + 2: dihedrals_p_count]

# Insert new dihedrals list

count = dihedrals_p_index + 2
for i in range(0, len(dihedrals)):
    b.insert(count, dihedrals[i] + "\n")
    count += 1

# [ dihedrals ] ; impropers

# In places where the double bond b/w c1 and c2 has changed to a single bond, the improper dihedrals need to be removed
# generate a list of impropers of interest (ones that will potentially be eliminated)

impropers_list = []
for i in range(0, tot_monomers):
    impropers_list.append([atoms*i + 26, atoms*i + 27, atoms*i + 86, atoms*i + 87])
    impropers_list.append([atoms*i + 41, atoms*i + 42, atoms*i + 111, atoms*i + 112])
    impropers_list.append([atoms*i + 56, atoms*i + 57, atoms*i + 136, atoms*i + 137])
    impropers_list.append([atoms*i + 25, atoms*i + 26, atoms*i + 27, atoms*i + 85])
    impropers_list.append([atoms*i + 40, atoms*i + 41, atoms*i + 42, atoms*i + 110])
    impropers_list.append([atoms*i + 55, atoms*i + 56, atoms*i + 57, atoms*i + 135])

imp_of_interest = []
for i in range(0, len(impropers_list)):
    if impropers_list[i][0] in c1 or impropers_list[i][1] in c1 or impropers_list[i][2] in c1 or impropers_list[i][3] in c1:
        # check against c1 since it is in both dihedrals of interest. c2 is involved in an improper that we want to keep
        imp_of_interest.append(impropers_list[i])

dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
while b[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
    dihedrals_imp_index += 1

dihedrals_imp_count = dihedrals_imp_index + 2
start_imp = dihedrals_imp_count

dihedrals_imp = []
while b[dihedrals_imp_count] != '\n':  # This is the last section in the input .itp file
    a_imp = [int(b[dihedrals_imp_count][0:6]), int(b[dihedrals_imp_count][6:13]), int(b[dihedrals_imp_count][13:20]), int(b[dihedrals_imp_count][20:27])]
    a_imp.sort()
    if a_imp not in imp_of_interest:
        dihedrals_imp.append(b[dihedrals_imp_count])
    dihedrals_imp_count += 1

# eliminate the old improper dihedrals list

del b[dihedrals_imp_index + 2: dihedrals_imp_count]

# Insert new improper dihedrals list

count = dihedrals_imp_index + 2
for i in range(0, len(dihedrals_imp)):
    b.insert(count, dihedrals_imp[i])
    count += 1

# [ virtual_sites4 ]

vsite_index = 0  # find index where [ virutal_sites4 ] section begins (propers)
while b[vsite_index].count('[ virtual_sites4 ]') == 0:
    vsite_index += 1

nv = 0
vsite_count = vsite_index + 2

for i in range(vsite_count, len(b)):  # This is the last section in the input .itp file
    vsite_count += 1
    nv += 1

# we need to make a new list of virtual sites that does not include in the sites which have been turned to real atoms
# term and H_new1 contain the indices of the hydrogens which need to be removed from the virtual sites list

virtual_sites = []
for i in range(vsite_index + 2, vsite_count):
    site = int(b[i][0:8])
    if site not in H_new2 and site not in H_new1:  # if the site is in the old virtual sites but not in H_new1 or tern
        virtual_sites.append(b[i])  # then add it to virtual_sites since those hydrogens are still dummies

# eliminate the old virtual_sites4 list

del b[vsite_index + 2: vsite_count]

# Insert new virtual sites list

count = vsite_index + 2
for i in range(0, len(virtual_sites)):
    b.insert(count, virtual_sites[i])
    count += 1

# Now write the new topology to a new file

target = open('crosslinked_new.itp', 'w')

for line in b:
    target.write(line)

target.close()

end = time.time()
print 'Total time elapsed: %s seconds' %(end - start)
print c1
print c2
print radical_c2
count = 0
for i in range(0, len(c2)):
    if c2[i] in radical_c2:
        count += 1

print 'Duplicates: %s' %count

print term
print term_c1
print other_c1
print H_new1
print H_new2