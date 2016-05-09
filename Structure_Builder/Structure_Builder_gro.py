# Script to build an array of LLC membrane membrane pores in .gro format
# It is necessary to use this script when there are more than 100,000 atoms

import numpy as np
import math

# Row at top of .gro file: (edit as necessary)
print 'This is a .gro file'

no_monomers = 6  # number of monomers packed per layer around a pore
no_atoms = 138  # number of atoms in one monomer excluding sodium ion
pore_radius = 3  # Radius of pore (Angstroms)
no_pores = 4  # number of pores to be simulated
dist_bw = 35  # distance between pores (Angstroms)
no_layers = 20  # Number of layers in a pore
sys_atoms = no_layers*no_monomers*no_pores*no_atoms  # total number of atoms in the system
print '%s' %sys_atoms
dist = 10  # distance between layers (Angstroms)
lines_of_text = 4  # lines of text at top of .pdb input file

f = open("Monomer_Configurations/monomer4.pdb", "r")

a = []
for line in f:
    a.append(line)

x_values_inp = []  # list to hold input values of x stored from .pdb file
y_values_inp = []  # list to hold input values of y stored from .pdb file
z_values_inp = []  # list to hold input values of z stored from .pdb file
positions_inp = []  # holds x, y, z coordinates of input .pdb file
identity = []  # holds the names of atom in the order that they appear in the .pdb file
# read specific entries in a text file
for i in range(lines_of_text, lines_of_text + no_atoms):  # searches relevant lines of text in file, f, being read
    # There are 137 atoms excluding sodium
    x_values_inp.append(float(a[i][26:38]))  # Use this to read specific entries in a text file
    y_values_inp.append(float(a[i][38:46]))  # makes sure I backtrack far enough to get all digits(i.e.38 instead of 42)
    z_values_inp.append(float(a[i][46:54]))
    positions_inp.append([x_values_inp[i - lines_of_text], y_values_inp[i - lines_of_text], z_values_inp[i - lines_of_text]])
    identity.append(a[i][13:16])  # hold name of atom (C, H or O)

# Get rid of whitespace in identities
for i in range(0, len(identity)):
    identity[i] = identity[i].strip()

# matrix to translate molecule to origin based on the position of atom 10 (Carbonyl carbon coming off benzene)
translation = np.matrix([[1, 0, 0,-(pore_radius + positions_inp[9][0])], [0, 1, 0,-(pore_radius+ positions_inp[9][1])],\
                         [0, 0, 1, -(pore_radius + positions_inp[9][2])], [0, 0, 0, 1]])
for i in range(0, len(positions_inp)):
    positions_inp[i].append(1)
    x = np.dot(translation, np.array(positions_inp[i]))
    positions_inp[i] = [x[0, 0], x[0, 1], x[0, 2]]

positions1 = []  # will hold the x,y,z coordinates of each atom of monomer 1
positions2 = []  # will hold the x,y,z coordinates of each atom of monomer 2
positions3 = []  # will hold the x,y,z coordinates of each atom of monomer 3
positions4 = []  # will hold the x,y,z coordinates of each atom of monomer 4
positions5 = []  # will hold the x,y,z coordinates of each atom of monomer 5
positions6 = []  # will hold the x,y,z coordinates of each atom of monomer 6
positions7 = []  # here for the case where 7 monomers pack
positions = [positions1, positions2, positions3, positions4, positions5, positions6, positions7]  # list of lists
x_values = []  # will hold x_values in the order that they appear in the positions matrix
y_values = []  # will hold y_values in the order that they appear in the positions matrix
z_values = []  # will hold z_values in the order that they appear in the positions matrix

# rotate coordinates and store each rotated coordinate as a separate list:

for k in range(0, len(positions_inp)):
    x = np.array(positions_inp[k])
    for i in range(1, no_monomers + 1):
        theta = i * math.pi / (no_monomers / 2.0)  # angle to rotate about axis determined from no of monomers per layer
        # Creates a rotation matrix to rotate input vector about y-axis making a new coordinate at evenly spaced points.
        # Each rotation belongs to a different monomer's position.The no of points corresponds to the number of monomers
        Rx = np.zeros((3, 3))  # makes a 3 x 3 zero matrix
        Rx[0, 0] = math.cos(theta)  # This line and subsequent edits to Rx fills in entries needed for rotation matrix
        Rx[1, 0] = math.sin(theta)
        Rx[0, 1] = -math.sin(theta)
        Rx[1, 1] = math.cos(theta)
        Rx[2, 2] = 1
        rot = np.dot(Rx, x)  # multiplies atomic coordinates by the rotation vector to generate new coordinates
        rot = [float(rot[0]), float(rot[1]), float(rot[2])]  # converts matrix to a list of floats
        positions[i - 1].append(rot)  # appends the atomic coordinates to 'positions'

for l in range(0, no_pores):  # loop to create multiple pores
    theta = 30  # angle which will be used to do hexagonal packing
    if l == 0:  # unmodified coordinates
        b = 0
        c = 0
    elif l == 1:  # move a pore directly down
        b = - 1
        c = 0
    elif l == 2:  # moves pore up and to the right
        b = math.cos(math.radians(90 - theta))
        c = -math.sin(math.radians(90 - theta))
    elif l == 3:  # moves a pore down and to the right
        b = -math.sin(math.radians(theta))
        c = -math.cos(math.radians(theta))
    for k in range(0, no_layers):
        for j in range(0, no_monomers):  # iterates over each monomer to create coordinates
            for i in range(0, len(positions[j]) - 1):  #
                x_values.append(b*dist_bw + positions[j][i][0])
                y_values.append(c*dist_bw + positions[j][i][1])
                z_values.append(k*dist + positions[j][i][2])
                if i + 1 + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l < 100000:
                    print '{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(1 + j + no_monomers*k + no_monomers*no_layers*l,
                        'LLC', identity[i], i + 1 + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l,
                        x_values[i+(no_atoms - 1)*j+(no_atoms - 1)*no_monomers*k+(no_atoms - 1)*no_monomers*no_layers*l]/10.0,
                        y_values[i + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l]/10.0,
                        z_values[i + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l]/10.0)
                else:
                    print '{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(1 + j + no_monomers*k + no_monomers*no_layers*l,
                        'LLC', identity[i], i + 1 + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l - 100000,
                        x_values[i+(no_atoms - 1)*j+(no_atoms - 1)*no_monomers*k+(no_atoms - 1)*no_monomers*no_layers*l]/10.0,
                        y_values[i + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l]/10.0,
                        z_values[i + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l]/10.0)


# Ions:

for l in range(0, no_pores):  # loop to create multiple pores
    theta = 30  # angle which will be used to do hexagonal packing
    if l == 0:  # unmodified coordinates
        b = 0
        c = 0
    elif l == 1:  # move a pore directly down
        b = - 1
        c = 0
    elif l == 2:  # moves pore up and to the right
        b = math.cos(math.radians(90 - theta))
        c = -math.sin(math.radians(90 - theta))
    elif l == 3:  # moves a pore down and to the right
        b = -math.sin(math.radians(theta))
        c = -math.cos(math.radians(theta))
    for k in range(0, no_layers):
        for j in range(0, no_monomers):  # iterates over each monomer to create coordinates
            x_values.append(b*dist_bw + positions[j][no_atoms - 1][0])
            y_values.append(c*dist_bw + positions[j][no_atoms - 1][1])
            z_values.append(k*dist + positions[j][no_atoms -1][2])
            if 1 + (no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers < 100000:
                print '{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.\
                    format(1 + j + no_monomers*k + no_monomers*no_layers*l + no_monomers*no_layers*no_pores,
                    identity[no_atoms - 1], identity[no_atoms - 1], 1 + (no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers,
                    x_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers]/10.0,
                    y_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers]/10.0,
                    z_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers]/10.0)
            else:
                print '{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.\
                    format(1 + j + no_monomers*k + no_monomers*no_layers*l + no_monomers*no_layers*no_pores,
                    identity[no_atoms - 1], identity[no_atoms - 1], 1 + (no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers - 100000,
                    x_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers]/10.0,
                    y_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers]/10.0,
                    z_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers]/10.0)

print '   0.00000   0.00000  0.00000'
