#!/usr/bin/python

# A script to clean up the topology file after completing crosslinking by getting rid of dummy atoms

import argparse

parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')  # allow input from user

parser.add_argument('-t', '--topology', default='crosslinked_new.itp', help = 'Name of input file')
parser.add_argument('-g', '--gro', default='wiggle.gro', help = '.gro file to be cleaned')
parser.add_argument('-m', '--monomer', default='HII', help = 'Residue name of monomer being used in the .gro')

args = parser.parse_args()

# Read the topology into a file

f = open(args.topology, 'r')

top = []
for line in f:
    top.append(line)

# Find the section containing all of the atoms

atoms_index = 0  # find index where [ atoms ] section begins
while top[atoms_index].count('[ atoms ]') == 0:
    atoms_index += 1

atoms_count = atoms_index + 2
nr = 0  # number of lines in 'atoms' section
while top[atoms_count] != '\n':
    atoms_count += 1  # increments the while loop
    nr += 1  # counts number of atoms

# Save only the atoms that are not dummy atoms

no_dummies = []  # new topology lines will be stored here
dummy_indices = []  # list to keep track of the dummy atoms being eliminated. This will be handy later
count = 1  # This will count only atoms that will be kept
count_all = 1  # This counts every atom that is being read
new_no = {}  # A dictionary to store new numbers as they are related to their old numbers (the dictionary keys are old numbers)
for i in range(atoms_index + 2, atoms_count):  # Renumber all of the atoms
    if top[i].count('hc_d') == 0:  # anything that isn't a dummy atom type will be kept and the number changed
        no_dummies.append(top[i][0:5].replace(top[i][0:5], '{:>5}'.format(str(count))) + top[i][5:30] + top[i][30:35].replace(top[i][30:35], '{:>5}'.format(str(count))) + top[i][35:len(top[i])])
        new_no[count_all] = count  # create dictionary entries
        count += 1
    else:
        dummy_indices.append(count_all)
    count_all += 1
# Now renumber all of the sections

# [ bonds ]

bonds_index = 0  # find index where [ bonds ] section begins
while top[bonds_index].count('[ bonds ]') == 0:
    bonds_index += 1

no_dummies.append('\n')  # add space between sections
for i in range(0, 2):  # add section headers and column labels
    no_dummies.append(top[bonds_index + i])

bond_count = bonds_index + 2
while top[bond_count] != '\n':  # replace all current numbers which are stored as keys with their dictionary value
    no_dummies.append(top[bond_count][0:6].replace(top[bond_count][0:6], '{:>6}'.format(str(new_no[int(top[bond_count][0:6])]))) + \
        top[bond_count][6:13].replace(top[bond_count][6:13], '{:>7}'.format(str(new_no[int(top[bond_count][6:13])]))) + top[bond_count][13:len(top[bond_count])])
    bond_count += 1  # increments while loop

# Find where the [ pairs ] section is located

pairs_index = 0  # find index where [ pairs ] section begins
while top[pairs_index].count('[ pairs ]') == 0:
    pairs_index += 1

no_dummies.append('\n')  # add space between sections
for i in range(0, 2):  # add section header and column labels
    no_dummies.append(top[pairs_index + i])

pairs_count = pairs_index + 2  # keep track of index of a
while top[pairs_count] != '\n':
    no_dummies.append(top[pairs_count][0:6].replace(top[pairs_count][0:6], '{:>6}'.format(str(new_no[int(top[pairs_count][0:6])]))) + \
        top[pairs_count][6:13].replace(top[pairs_count][6:13], '{:>7}'.format(str(new_no[int(top[pairs_count][6:13])]))) + top[pairs_count][13:len(top[pairs_count])])
    pairs_count += 1

# [ angles ]

angles_index = 0  # find index where [ angles ] section begins
while top[angles_index].count('[ angles ]') == 0:
    angles_index += 1

no_dummies.append('\n')  # add space between sections
for i in range(0, 2):  # add section header and column labels
    no_dummies.append(top[angles_index + i])

angle_count = angles_index + 2  # keep track of index of a
while top[angle_count] != '\n':
    no_dummies.append(top[angle_count][0:6].replace(top[angle_count][0:6], '{:>6}'.format(str(new_no[int(top[angle_count][0:6])]))) + \
        top[angle_count][6:13].replace(top[angle_count][6:13], '{:>7}'.format(str(new_no[int(top[angle_count][6:13])]))) + \
        top[angle_count][13:20].replace(top[angle_count][13:20], '{:>7}'.format(str(new_no[int(top[angle_count][13:20])]))) + \
        top[angle_count][20:len(top[angle_count])])
    angle_count += 1

# [ dihedrals ] ; propers

dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
while top[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
    dihedrals_p_index += 1

no_dummies.append('\n')  # add space between sections
for i in range(0, 2):  # add section header and column labels
    no_dummies.append(top[dihedrals_p_index + i])

dihedrals_p_count = dihedrals_p_index + 2  # keep track of index of a
while top[dihedrals_p_count] != '\n':
    no_dummies.append(top[dihedrals_p_count][0:6].replace(top[dihedrals_p_count][0:6], '{:>6}'.format(str(new_no[int(top[dihedrals_p_count][0:6])]))) + \
        top[dihedrals_p_count][6:13].replace(top[dihedrals_p_count][6:13], '{:>7}'.format(str(new_no[int(top[dihedrals_p_count][6:13])]))) + \
        top[dihedrals_p_count][13:20].replace(top[dihedrals_p_count][13:20], '{:>7}'.format(str(new_no[int(top[dihedrals_p_count][13:20])]))) + \
        top[dihedrals_p_count][20:27].replace(top[dihedrals_p_count][20:27], '{:>7}'.format(str(new_no[int(top[dihedrals_p_count][20:27])]))) + \
        top[dihedrals_p_count][27:len(top[dihedrals_p_count])])
    dihedrals_p_count += 1

# [ dihedrals ] ; impropers

dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
while top[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
    dihedrals_imp_index += 1

no_dummies.append('\n')  # add space between sections
for i in range(0, 2):  # add section header and column labels
    no_dummies.append(top[dihedrals_imp_index + i])

dihedrals_imp_count = dihedrals_imp_index + 2
while top[dihedrals_imp_count] != '\n':
    no_dummies.append(top[dihedrals_imp_count][0:6].replace(top[dihedrals_imp_count][0:6], '{:>6}'.format(str(new_no[int(top[dihedrals_imp_count][0:6])]))) + \
        top[dihedrals_imp_count][6:13].replace(top[dihedrals_imp_count][6:13], '{:>7}'.format(str(new_no[int(top[dihedrals_imp_count][6:13])]))) + \
        top[dihedrals_imp_count][13:20].replace(top[dihedrals_imp_count][13:20], '{:>7}'.format(str(new_no[int(top[dihedrals_imp_count][13:20])]))) + \
        top[dihedrals_imp_count][20:27].replace(top[dihedrals_imp_count][20:27], '{:>7}'.format(str(new_no[int(top[dihedrals_imp_count][20:27])]))) + \
        top[dihedrals_imp_count][27:len(top[dihedrals_imp_count])])
    dihedrals_imp_count += 1

f = open('crosslinked.itp', 'w')
for i in range(0, atoms_index + 2):
    f.writelines(top[i])
for i in no_dummies:
    f.writelines([i])

# Now we need to rewrite the .gro file

f = open(args.gro, 'r')
gro = []
for line in f:
    gro.append(line)

# Find the beginning of the file
start = 0
while gro[start].count(args.monomer) == 0:
    start += 1

no_dummy_gro = []

for i in range(start, start + count_all - 1):  # looks through all indices of monomer
    if int(gro[i][15:20]) not in dummy_indices:
        no_dummy_gro.append(gro[i][0:15] + gro[i][15:20].replace(gro[i][15:20], '{:>5}'.format(str(new_no[int(gro[i][15:20])]))) + gro[i][20:len(gro[i])])

other_index = count_all - 1 + start
new_index_count = len(no_dummy_gro) + 1
for i in range(other_index, len(gro) - 1):
    no_dummy_gro.append(gro[i][0:15] + gro[i][15:20].replace(gro[i][15:20], '{:>5}'.format(str(new_index_count)) + gro[i][20:len(gro[i])]))
    other_index += 1
    new_index_count += 1

f = open('wiggle_no_dummies.gro', 'w')
f.writelines(['%s' %gro[0], '%s' %len(no_dummy_gro) + '\n'])
for i in no_dummy_gro:
    f.writelines(i)
f.writelines(['%s' %gro[len(gro) - 1]])