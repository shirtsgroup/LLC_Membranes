#!/usr/bin/python

# Characterize Crosslinked systems

import argparse
import math

parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')  # allow input from user

# Flags
parser.add_argument('-g', '--gro', default='wiggle_no_dummies.gro', help = 'Name of the final, cleaned .gro file after \
                                                                            crosslinking finishes')
parser.add_argument('-t', '--top', default='crosslinked.itp', help = '.itp file associated with input .gro file')
parser.add_argument('-m', '--monomer', default='HII', help = 'Type of monomer used in crosslinked system')
parser.add_argument('-p', '--no_pores', default=4, help = 'Number of Pores')
args = parser.parse_args()

f = open(args.gro, 'r')

g = []
for line in f:
    g.append(line)

f = open(args.top, 'r')

t = []
for line in f:
    t.append(line)

exec 'from Scripts.LC_class import %s' % args.monomer
exec 'c1_atoms = %s.c1_atoms' % args.monomer
exec 'c2_atoms = %s.c2_atoms' % args.monomer
exec 'tails = %s.tails' % args.monomer

# Find which atoms the carbons are bonded to

# find [ atoms ] section in topology

atoms_index = 0
while t[atoms_index].count('[ atoms ]') == 0:
    atoms_index += 1

count = atoms_index + 2
while t[count] != '\n':
    count += 1

c1_no = []
c2_no = []

for i in range(atoms_index + 2, count):
    if str.strip(t[i][23:28]) in c1_atoms:
        c1_no.append(int(t[i][0:5]))
    if str.strip(t[i][23:28]) in c2_atoms:
        c2_no.append(int(t[i][0:5]))

# Create a dictionary with keys that are atom numbers and values which are atom types. This will be useful later
atom_types = {}
for i in range(atoms_index + 2, count):
    atom_types[str.strip(t[i][0:5])] = str.strip(t[i][23:28])

# find the [ bonds ] section

bonds_index = 0
while t[bonds_index].count('[ bonds ]') == 0:
    bonds_index += 1

bond_count = bonds_index + 2
while t[bond_count] != '\n':
    bond_count += 1

c1_bonds = []
c1_xlinks = []
for i in range(0, len(c1_no)):
    c1_bonds.append([])
    c1_xlinks.append([])

for i in range(bonds_index + 2, bond_count):
    a = int(t[i][0:6])
    b = int(t[i][6:13])
    if a in c1_no:
        c1_bonds[c1_no.index(a)].append(b)
    if b in c1_no:
        c1_bonds[c1_no.index(b)].append(a)

# Filter out regular bonds -- i.e. bonds that are not crosslinked but are just a part of the chain
# Do this by checking if the atom number is a c1 or c2

for i in range(0, len(c1_bonds)):
    for j in range(0, len(c1_bonds[i])):
        if atom_types[str(c1_bonds[i][j])] in c2_atoms:
            if c1_bonds[i][j] != c2_no[i]:
                c1_xlinks[i].append(c1_bonds[i][j])

# Find the total number of crosslinks -- if a list entry in c1_xlinks is empty then there is no crosslink there

xlinks = 0
for i in range(0, len(c1_xlinks)):
    if c1_xlinks[i]:
        xlinks += 1

# Let's find the range of atom numbers for each monomer
# Each atom in the monomer has a unique atom name (C, C1, C2, C3,...,O, 01 etc.). A new monomer can be found when one of
# atom names repeat themselves. First, find the very first atom name:

leading_atom = str.strip(t[atoms_index + 2][23:28])

mon_ranges = [[int(t[atoms_index + 2][0:5])]]
index = 0
for i in range(atoms_index + 3, count):
    if str.strip(t[i][23:28]) == leading_atom:
        mon_ranges.append([])
        index += 1
    mon_ranges[index].append(int(t[i][0:5]))

# Now find the number of crosslinks that are within the same monomer by seeing if the crosslinked atom is in the same
# range as the one it is bonded to

intra_mon = 0
for i in range(0, len(c1_xlinks)):
    for j in range(0, len(c1_xlinks[i])):
        if c1_xlinks[i][j] in mon_ranges[int(math.floor(i/tails))]:
            intra_mon += 1

# Find the distance between crosslinked atoms

distances = []

top_lines = 0
while g[top_lines].count('%s' %args.monomer) == 0:
    top_lines += 1

for i in range(0, len(c1_xlinks)):
    if c1_xlinks[i]:
        c1_x = float(g[c1_no[i] + top_lines - 1][20:28])
        c1_y = float(g[c1_no[i] + top_lines - 1][28:36])
        c1_z = float(g[c1_no[i] + top_lines - 1][36:44])
        c2_x = float(g[c1_xlinks[i][0] + top_lines - 1][20:28])
        c2_y = float(g[c1_xlinks[i][0] + top_lines - 1][28:36])
        c2_z = float(g[c1_xlinks[i][0] + top_lines - 1][36:44])
        distances.append(math.sqrt((c1_x - c2_x)**2 + (c1_y - c2_y)**2 + (c1_z - c2_z)**2))

bonds_across_pbcs = 0
for i in range(0, len(distances)):
    if distances[i] > 0.3:
        bonds_across_pbcs += 1

# Calculate crosslinks between pores

pore_ranges = []  # a range of indices for each pore
mon_per_pore = len(c1_xlinks)/tails/args.no_pores

for i in range(0, args.no_pores):
    pore_ranges.append([])
    for j in range(0, mon_per_pore):
        for k in range(0, len(mon_ranges[i*mon_per_pore + j])):
            pore_ranges[i].append(mon_ranges[i*mon_per_pore + j][k])

inter_pore = 0
for i in range(0, len(c1_xlinks)):
    for j in range(0, len(c1_xlinks[i])):
        if c1_xlinks[i][j] not in pore_ranges[int(math.floor(i/tails/mon_per_pore))]:
            inter_pore += 1

# find crosslinks within same pore but not intra-monomer

intra_pore = 0

for i in range(0, len(c1_xlinks)):
    for j in range(0, len(c1_xlinks[i])):
        if c1_xlinks[i][j] in pore_ranges[int(math.floor(i/tails/mon_per_pore))]:
            if c1_xlinks[i][j] not in mon_ranges[int(math.floor(i/tails))]:
                intra_pore += 1

# Outputs

print 'Total Number of crosslinks: %s' %xlinks
print ''
print ' ______________________________________'
print '|____________________| Number |   %    |'
print '|Intra-Monomer xlinks|{:^8}|{:^8}|'.format(intra_mon, 100*intra_mon/xlinks)
# print '|Xlinks across PBCs  |{:^8}|{:^8}|'.format(bonds_across_pbcs, 100*bonds_across_pbcs/xlinks)
print '| Inter-Pore xlinks  |{:^8}|{:^8}|'.format(inter_pore, 100*inter_pore/xlinks)
print '| Intra-Pore xlinks  |{:^8}|{:^8}|'.format(intra_pore, 100*intra_pore/xlinks)
print ' --------------------------------------'