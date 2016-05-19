# Converts a gro file to a .xyz file

f = open("perfect.gro", "r")  # .gro file whose atomic coordinates will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

top_text = 2  # Lines at top of file which will be ignored

no_atoms = len(a) - top_text - 1  # every line of the .gro file has atomic coordinates except the top and bottom

print no_atoms
print '.gro file now in .xyz format'

x = []  # list to hold input values of x stored from .gro file
y = []  # list to hold input values of y stored from .gro file
z = []  # list to hold input values of z stored from .gro file
identity = []  # holds the names of atom in the order that they appear in the .pdb file
# read specific entries in a text file
for i in range(top_text, top_text + no_atoms):  # searches relevant lines of text in file, f, being read
    # There are 137 atoms excluding sodium
    x.append(float(a[i][21:30]))  # Use this to read specific entries in a text file
    y.append(float(a[i][30:38]))  # makes sure I backtrack far enough to get all digits(i.e.38 instead of 42)
    z.append(float(a[i][38:45]))
    identity.append(str(a[i][9:15]))  # hold name of atom (C, H or O)

for i in range(0, len(identity)):
    if identity[i].count('C') == 1:
        identity[i] = 'C'
    elif identity[i].count('H') == 1:
        identity[i] = 'H'
    elif identity[i].count('O') == 1:
        identity[i] = 'O'
    elif identity[i].count('NA') == 1:
        identity[i] = 'NA'
    else:
        print 'unrecognized atoms(s)'

for i in range(0, len(identity)):
    print '{:>2s}{:12}{:12}{:12}'.format(identity[i], x[i], y[i], z[i])

