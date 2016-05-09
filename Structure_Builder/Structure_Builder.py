# Script to build an array of LLC membrane membrane pores!

import numpy as np
import math
import datetime

# Rows at top of .pdb file: (edit as necessary)
Header = 'LLC_membrane'
Title = 'Single Pore'
Compound = 'Na_GA3C11'
Author = 'Ben Coscia'
Today = datetime.date.today()
# Format top of .pdb file:
# print '{:10s} {:59s} {:4s} {:4d}'.format('HEADER', Header, 'NONE', 1)
# print '{:10s} {:59s} {:4s} {:4d}'.format('TITLE', Title, 'NONE', 2)
# print '{:10s} {:59s} {:4s} {:4d}'.format('COMPND', Compound, 'NONE', 3)
# print '{:10s} {:59s} {:4s} {:4d}'.format('AUTHOR', Author, 'NONE', 4)
# print '{:10s} {:59s} {:4s} {:4d}'.format('REVDAT', str(Today), 'NONE', 5)

no_monomers = 6  # number of monomers packed per layer around a pore
no_atoms = 138  # number of atoms in one monomer excluding sodium ion
pore_radius = 3  # Radius of pore (unsure of units right now)
no_pores = 4  # number of pores to be simulated
dist_bw = 40  # distance between pores (units tbd)
no_layers = 20  # Number of layers in a pore
dist = 10 # distance between layers (units tbd)
lines_of_text = 5  # lines of text at top of .pdb file

f = open("monomer.pdb", "r")
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
#identity.append('NA')

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
                print '{:5s}{:6d}  {:4s} {:8d}   {:8.3f} {:8.3f}{:8.3f}  {:4.2f}  {:4.2f}          {:>2}{:2s}'\
                    .format('ATOM', i + 1 + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l,
                    identity[i], 0, x_values[i+(no_atoms - 1)*j+(no_atoms - 1)*no_monomers*k+(no_atoms - 1)*no_monomers*no_layers*l],
                    y_values[i + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l],
                    z_values[i + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l], 0.00,
                    0.00, 'C', '+0')

# Next, we need to extract connectivity information from the .pdb file
# this was a pain so appreciate it!

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
            print '{:5s}{:6d}  {:4s} {:8d}   {:8.3f} {:8.3f}{:8.3f}  {:4.2f}  {:4.2f}          {:>2}{:2s}'\
                .format('ATOM', 1 + (no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers,
                identity[no_atoms - 1], 0, x_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers],
                y_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers],
                z_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers], 0.00,
                0.00, 'C', '+0')

c1 = []  # connectivity column 1 from .pdb file
c2 = []  # connectivity column 2 from .pdb file
c3 = []  # connectivity column 3 from .pdb file
c4 = []  # connectivity column 4 from .pdb file
c5 = []  # connectivity column 5 from .pdb file

for l in range(0, no_pores):
    for k in range(0, no_layers):
        for j in range(0, no_monomers):
            for i in range(143, 203):  # loops through rows which contain connectivity information
                if int(a[i][6:11]) != 0:  # need to retain 0's since they mean that there is no connection at that point
                    c1.append(j*(no_atoms - 1) + int(a[i][6:11]))  # searches for up to a 5 digit number
                else:
                    c1.append(0)
                if int(a[i][11:16]) != 0:
                    c2.append(j*(no_atoms - 1) + int(a[i][11:16]))
                else:
                    c2.append(0)
                if int(a[i][16:21]) != 0:
                    c3.append(j*(no_atoms - 1) + int(a[i][16:21]))
                else:
                    c3.append(0)
                if int(a[i][21:26]) != 0:
                    c4.append(j*(no_atoms - 1) + int(a[i][21:26]))
                else:
                    c4.append(0)
                if int(a[i][26:31]) != 0:
                    c5.append(j*(no_atoms - 1) + int(a[i][26:31]))
                else:
                    c5.append(0)
                print '{:6s}{:5d}{:5d}{:5d}{:5d}{:5d}                                         {:4s}{:4d}'\
                    .format('CONECT', c1[j * 60 + i - 143], c2[j * 60 + i - 143], c3[j * 60 + i - 143],
                    c4[j * 60 + i - 143], c5[j * 60 + i - 143], 'NONE', j*60 + 5 + (no_atoms -1)*(no_monomers-1) + i + 1)
# line numbering (last column) is a bit off but it doesn't matter
print 'END'
