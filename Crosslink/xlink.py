# Script to crosslink LLC monomers based on distance between carbons at the end of a simulation

#!/usr/bin/python

import os
import argparse
import numpy as np
import math
import time
import matplotlib.pyplot as plt
start = time.time()
location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')

# Flags
parser.add_argument('-i', '--input', default='wiggle.gro', help = 'Name of input file')
parser.add_argument('-t', '--type', default='LLC', help='Type of monomer')
parser.add_argument('-l', '--layers', default=20, help='Number of layers')
parser.add_argument('-p', '--pores', default=4, help='Number of pores')
parser.add_argument('-n', '--no_monomers', default=6, help='Number of monomers per layer')
parser.add_argument('-c', '--cutoff', default=5, help='Cutoff distance for cross-linking. Bottom x % of the distribution'
                                                      'of distances will be cross-linked ')

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
    atoms = 137  # not including sodium
    atoms_list = ['O5', 'O6', 'O7', '08', 'O9', 'C32', 'C18', 'C46', 'C33', 'C19', 'C47', 'C34', 'C20', 'C47', 'H49',
                  'H24', 'H74', 'H50', 'H25', 'H75', 'H51', 'H26', 'H76']
    topology = 'HII_mon.itp'
    xlink_atoms = 6  # number of atoms involved in cross-linking

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

for i in range(0, len(C1x)):
    for k in range(0, len(C2x)):
        if exclude[k, i] == 1:  # make sure that C1 and C2 that are a part of the same monomer do not factor into this calculation
            dist[k, i] = 1000  # artificially high number to keep it from interfering
        else:
            dist[k, i] = math.sqrt((C1x[i] - C2x[k])**2 + (C1y[i] - C2y[k])**2 + (C1z[i] - C2z[k])**2)
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

print len(min_list)
cutoff = min(min_list)  # The minimum value left after the modification of min_list is the cutoff value
print cutoff
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
end = time.time()
print 'Total time elapsed: %s seconds' %(end - start)

# reduce lists down to unique atom pairs


def uniq(input1, input2):
    output1 = []
    output2 = []
    for i in range(0, len(input1)):
        if input1[i] not in output1 and input2[i] not in output2:
            output1.append(input1[i])
            output2.append(input2[i])
    return output1, output2

c1, c2 = uniq(C1_no, C2_no)
print len(c1)
print len(c2)

dist_list = []
for i in range(0, len(c1)):
    C1 = int(c1[i]) + 1
    C2 = int(c2[i]) + 1
    print C1
    print C2
    dist_list.append(math.sqrt((float(a[C1][20:28]) - float(a[C2][20:28]))**2 +
                     (float(a[C1][28:36]) - float(a[C2][28:36]))**2 +
                     (float(a[C1][36:44]) - float(a[C2][36:44]))**2))

for i in range(0, len(dist_list)):
    if dist_list[i] > 0.4:
        print '%s), %s, %s, %s' %(i, c1[i], c2[i], dist_list[i])