# Delete parameters for bonds, angles and dihedrals that acpype.py fills in automatically

f = open("HII_mon.itp", "r")

a = []
for line in f:
    a.append(line)

angles_index = 0  # find index where [ angles ] section begins
while a[angles_index].count('[ angles ]') == 0:
    angles_index += 1

na = 0  # number of lines in the 'angles' section
angle_count = angles_index + 2  # keep track of index of a
while a[angle_count] != '\n':
    angle_count += 1
    na += 1

dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
while a[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
    dihedrals_p_index += 1

ndp = 0  # number of lines in the 'dihedrals ; proper' section
dihedrals_p_count = dihedrals_p_index + 3  # keep track of index of a
while a[dihedrals_p_count] != '\n':
    dihedrals_p_count += 1
    ndp += 1

dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
while a[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
    dihedrals_imp_index += 1

ndimp = 0  # number of lines in the 'dihedrals ; impropers' section
dihedrals_imp_count = dihedrals_imp_index + 3

for i in range(dihedrals_imp_count, len(a)):  # This is the last section in the input .itp file
    dihedrals_imp_count += 1
    ndimp += 1

bonds_index = 0  # find index where [ bonds ] section begins
while a[bonds_index].count('[ bonds ]') == 0:
    bonds_index += 1

nb = 0  # number of lines in the 'bonds' section
bond_count = bonds_index + 2
while a[bond_count] != '\n':
    bond_count += 1  # increments while loop
    nb += 1  # counting number of lines in 'bonds' section

for i in range(0, bonds_index + 2):
    print a[i],

for i in range(bonds_index + 2, bond_count):
    print a[i][0:17]

for i in range(bond_count, angles_index + 2):
    print a[i],

for i in range(angles_index + 2, angle_count):
    print a[i][0:28]

for i in range(angle_count, angle_count + 4):
    print a[i],

for i in range(dihedrals_p_index + 3, dihedrals_p_count):
    print a[i][0:35]

for i in range(dihedrals_p_count, dihedrals_p_count + 4):
    print a[i],

for i in range(dihedrals_imp_index + 3, dihedrals_imp_count):
    print a[i][0:35]


