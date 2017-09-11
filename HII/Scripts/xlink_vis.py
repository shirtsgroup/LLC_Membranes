#!/usr/bin/python

"""
Create a .pdb file which contains coordinates of crosslinked atoms. To make it less busy, each crosslink is represented
by a coordinate halfway between the bonded carbons.
"""

import argparse
from Get_Positions import get_positions
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Create .gro for crosslink visualization')

parser.add_argument('-f', '--file', default='wiggle_no_dummies.gro', help='Coordinate file in .gro format')
parser.add_argument('-t', '--top', default='crosslinked.itp', help='Final crosslinked topology')
parser.add_argument('-c', '--c1', default=['C20', 'C34', 'C48'])
parser.add_argument('-C', '--c2', default=['C19', 'C33', 'C47'])
parser.add_argument('-s', '--subtitle', default='')

args = parser.parse_args()

pos, _, _, box = get_positions('%s' % args.file, 'sys', 'HII', 'no')  # only get the first return values which is a position array

f = open('%s' % args.top, 'r')

a = []
for line in f:
    a.append(line)
f.close()

# find the [ atoms ] section
atoms_index = 0  # find index where [ atoms ] section begins
while a[atoms_index].count('[ atoms ]') == 0:
    atoms_index += 1

# find the [ bonds ] section
bonds_index = 0  # find index where [ bonds ] section begins
while a[bonds_index].count('[ bonds ]') == 0:
    bonds_index += 1
# find the end of the [ bonds ] section
end_bonds = bonds_index
while a[end_bonds] != '\n':
    end_bonds += 1

# find the crosslinked and/or terminated carbons
c1 = []
c2 = []

for i in range(atoms_index, bonds_index):
    if str.strip(a[i][23:28]) in args.c1 and a[i][-2] == 'T':  # I must've used a space at the end so its a[i][-2] instead of [-1]
        c1.append(int(a[i][0:5]))
    elif str.strip(a[i][23:28]) in args.c2 and a[i][-2] == 'T':
        c2.append(int(a[i][0:5]))

# Now look through the bonds section to see between which carbons the crosslinks occured
c1x = []
c2x = []

for i in range(bonds_index + 2, end_bonds):
    b1 = int(a[i][0:6])
    b2 = int(a[i][6:13])
    if b1 in c1:
        if b2 in c2 and b2 != b1 + 1:  # make sure that we aren't counting adjacent carbons in the same monomer
            c1x.append(b1)
            c2x.append(b2)
    if b2 in c1:  # This case actually doesn't occur based on how xlink.py rewrites the topology
        if b1 in c2 and b1 != b2 - 1:
            c1x.append(b2)
            c2x.append(b1)

# Get positions for each of crosslinking carbon
xlinks = len(c1x)

c1_coords = np.zeros([3, xlinks])
c2_coords = np.zeros([3, xlinks])

for i in range(xlinks):
    c1_coords[:, i] = pos[:, c1x[i] - 1, 0]
    c2_coords[:, i] = pos[:, c2x[i] - 1, 0]

# Now create artificial coordinates for a new .gro file highlighting the crosslinks by averaging the locations of each
# bonded carbon


art_coords = np.zeros([3, xlinks])

for i in range(xlinks):
    art_coords[:, i] = (c1_coords[:, i] + c2_coords[:, i])/2.

box_f = '{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}'.format(float(box[0]), float(box[1]), float(box[2]),
                                                                        float(box[3]), float(box[4]), float(box[5]),
                                                                        float(box[6]), float(box[7]), float(box[8]))
f = open('xlinks.gro', 'w')

f.write('This is a .gro containing average coordinates of crosslinked carbons\n')
f.write('%s\n' % xlinks)
for i in range(xlinks):
    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(i + 1, 'NA', 'NA', i + 1, art_coords[0, i],
                                                                                    art_coords[1, i], art_coords[2, i]))
f.write(box_f)
f.close()

# find radial distribution of ions

na_pos = get_positions('%s' % args.file, ['NA'], 'HII', 'no')[0]
no_pores = 4
no_ions = len(na_pos[0, :, 0])
ion_ppore = no_ions/no_pores

pores = np.zeros([no_pores, 3, ion_ppore])  # [pore numbers, xyz dimensions, ion number within pore]

for p in range(no_pores):
    for i in range(ion_ppore*p, ion_ppore*(p + 1)):
        pores[p, :, i - p*ion_ppore] = na_pos[:, i, 0]

# find the average x-y positions of the pore centers
centers = np.zeros([no_pores, 2])
for p in range(no_pores):
    for d in range(2):
        centers[p, d] = np.mean(pores[p, d, :])

# Use the pore center as a central axis and create a radial distribution from that
rdf = np.zeros([no_ions])
for p in range(no_pores):
    for i in range(ion_ppore):
        dist = np.linalg.norm(pores[p, 0:2, i]-centers[p, :])
        rdf[p*ion_ppore + i] = dist

print np.mean(rdf)

bins = np.linspace(0, 1.5, 50)
plt.hist(rdf, bins=bins)
if len(args.subtitle) > 0:
    plt.suptitle('Radial Distribution of Sodium Ions')
    plt.title('%s' % args.subtitle)
else:
    plt.title('Radial Distribution of Sodium Ions')
plt.xlabel('Radial distance from pore center (nm)')
plt.ylabel('Count')
plt.show()
