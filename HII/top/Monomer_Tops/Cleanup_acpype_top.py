#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser(description = 'Build LLC Structure')

parser.add_argument('-f', '--file', help = 'File to be cleaned')

args = parser.parse_args()

f = open('%s.itp' %args.file, 'r')

b = []
for line in f:
    b.append(line)

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
    bond_count += 1  # increments while loop

# [ pairs ]

pairs_index = 0
while b[pairs_index].count('[ pairs ]') == 0:
	pairs_index += 1

pairs_count = pairs_index + 2
while b[pairs_count] != '\n':
	pairs_count += 1

# [ angles ]

angles_index = 0  # find index where [ angles ] section begins
while b[angles_index].count('[ angles ]') == 0:
    angles_index += 1

na = 0  # number of lines in the 'angles' section
angle_count = angles_index + 2  # keep track of index of a
angles_prev = []
while b[angle_count] != '\n':
    angle_count += 1

# [ dihedrals ] ; propers

dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
while b[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
    dihedrals_p_index += 1

ndp = 0  # number of lines in the 'dihedrals ; proper' section
dihedrals_p_count = dihedrals_p_index + 2  # keep track of index of a
dihedrals_prev = []
while b[dihedrals_p_count] != '\n':
    dihedrals_p_count += 1

dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
while b[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
    dihedrals_imp_index += 1

ndimp = 0  # number of lines in the 'dihedrals ; impropers' section
dihedrals_imp_count = dihedrals_imp_index + 3
start_imp = dihedrals_imp_count

dihedrals_imp = []
#for i in range(dihedrals_imp_count, len(b)):
while b[dihedrals_imp_count] != '\n':  # This is the last section in the input .itp file
    dihedrals_imp_count += 1

a = []
for i in range(0, atoms_index + 1):
	a.append(str.strip(b[i]))

for i in range(atoms_index + 1, atoms_count):
	a.append(b[i][0:60])

for i in range(atoms_count, bonds_index + 1):
	a.append(str.strip(b[i]))

for i in range(bonds_index + 1, bond_count - 6):
	a.append(b[i][0:20])

for i in range(bond_count -6, bond_count):
	a.append('    '+ str.strip(b[i][0:20]))

for i in range(bond_count, pairs_index + 2):
	a.append(str.strip(b[i]))

for i in range(pairs_index + 2, pairs_count):
	a.append(b[i][0:20])

for i in range(pairs_count, angles_index + 1):
	a.append(str.strip(b[i]))

for i in range(angles_index + 1, angle_count):
	a.append(b[i][0:29])

for i in range(angle_count, dihedrals_p_index + 1):
	a.append(str.strip(b[i]))

for i in range(dihedrals_p_index + 2, dihedrals_p_count):
	a.append(b[i][0:35])

for i in range(dihedrals_p_count, dihedrals_imp_index + 1):
	a.append(str.strip(b[i]))

for i in range(dihedrals_imp_index + 2, dihedrals_imp_count):
	a.append(b[i][0:35])

for i in range(dihedrals_imp_count, len(b)):
	a.append(str.strip(b[i]))

#for i in range(0, len(a)):
#	print a[i]

f = open('%s_cleaned.itp' %args.file, 'w')

for line in a:
	f.writelines(line + '\n')
