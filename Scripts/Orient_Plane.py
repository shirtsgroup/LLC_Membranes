#!/usr/bin/python

# Orient plane with origin
import numpy as np
import math
import os
import argparse


location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

parser = argparse.ArgumentParser(description='Build LLC Structure')
parser.add_argument('-t', '--type', default='HII', type = str, help = 'membrane type')
parser.add_argument('-o', '--out', default='initial.gro', help='name of output file')
parser.add_argument('-i', '--input', default='monomer2.pdb', help = 'Path to input file')
parser.add_argument('-l', '--layers', default=20, type=int, help = 'Number of Layers')
parser.add_argument('-m', '--monomers', default=6, type=int, help = 'Monomers per layer')
parser.add_argument('-r', '--radius', default=6, type=float, help = 'Initial Pore Radius (Angstroms)')
parser.add_argument('-p', '--p2p', default=45, type=float, help = 'Initial Pore to Pore Distance')
parser.add_argument('-n', '--nopores', default=4, type=int, help = 'Number of Pores')
parser.add_argument('-d', '--dbwl', default=5, type=float, help = 'Distance between layers')
parser.add_argument('-s', '--layer_distribution', default='uniform', help = 'The distribution of monomes per layer')
parser.add_argument('-a', '--alt_1', default=6, type=int, help='Monomers per layer for the first type of alternating layer')
parser.add_argument('-A', '--alt_2', default=8, type=int, help='Monomers per layer for the second type of alternating layer')
args = parser.parse_args()


def read_pdb_coords(file):

    a = []
    for line in file:
        a.append(line)
    f.close()

    no_atoms = 0  # number of atoms in one monomer including sodium ion
    for i in range(0, len(a)):
        no_atoms += a[i].count('ATOM')

    lines_of_text = 0  # lines of text at top of .pdb input file
    for i in range(0, len(a)):
        if a[i].count('ATOM') == 0:
            lines_of_text += 1
        if a[i].count('ATOM') == 1:
            break

    xyz = np.zeros([3, no_atoms])
    identity = np.zeros([no_atoms], dtype=object)
    for i in range(lines_of_text, lines_of_text + no_atoms):  # searches relevant lines of text in file, f, being read
        # There are 137 atoms excluding sodium
        xyz[:, i - lines_of_text] = [float(a[i][26:38]), float(a[i][38:46]), float(a[i][46:54])]  # Use this to read specific entries in a text file
        identity[i - lines_of_text] = str.strip(a[i][12:16])

    return xyz, identity, no_atoms, lines_of_text


def read_gro_coords(file):

    a = []
    for line in file:
        a.append(line)
    f.close()

    lines_of_text = 2
    no_atoms = len(a) - lines_of_text - 1  # subtract one for the bottom box vector line

    xyz = np.zeros([3, no_atoms])
    identity = np.zeros([no_atoms], dtype=object)
    for i in range(lines_of_text, lines_of_text + no_atoms):  # searches relevant lines of text in file, f, being read
        xyz[:, i - lines_of_text] = [float(a[i][20:28])*10, float(a[i][28:36])*10, float(a[i][36:44])*10]
        identity[i - lines_of_text] = str.strip(a[i][11:16])

    return xyz, identity, no_atoms, lines_of_text

f = open("%s/../Structure-Files/HII_Monomer_Configurations/%s" % (location, args.input), "r")
t = 'HII'
no_ions = 1

if args.input.endswith('.pdb'):
    xyz, identity, no_atoms, lines_of_text = read_pdb_coords(f)
elif args.input.endswith('.gro'):
    xyz, identity, no_atoms, lines_of_text = read_gro_coords(f)
else:
    print 'Please input a valid file type (.gro or .pdb)'

layer_distribution = np.zeros([args.layers*args.nopores], dtype=int)

if args.layer_distribution == 'uniform':
    for i in range(0, len(layer_distribution)):
        layer_distribution[i] = args.monomers
if args.layer_distribution == 'alternating':
    for i in range(layer_distribution.shape[-1]):
        if i % 2 == 0:
            layer_distribution[i] = args.alt_1
        if i % 2 == 1:
            layer_distribution[i] = args.alt_2

no_monomers = args.monomers  # number of monomers packed per layer around a pore
pore_radius = args.radius  # Radius of pore (unsure of units right now)
no_pores = args.nopores  # number of pores to be simulated
dist_bw = args.p2p  # distance between pores (units tbd)
no_layers = args.layers  # Number of layers in a pore
sys_atoms = sum(layer_distribution)*no_atoms  # total number of atoms in the system
dist = args.dbwl  # distance between layers (units tbd)

# x_values_inp = []  # list to hold input values of x stored from .pdb file
# y_values_inp = []  # list to hold input values of y stored from .pdb file
# z_values_inp = []  # list to hold input values of z stored from .pdb file
# positions_inp = []  # holds x, y, z coordinates of input .pdb file
# identity = []  # holds the names of atom in the order that they appear in the .pdb file
# read specific entries in a text file
# for i in range(lines_of_text, lines_of_text + no_atoms):  # searches relevant lines of text in file, f, being read
#     # There are 137 atoms excluding sodium
#     x_values_inp.append(float(a[i][26:38]))  # Use this to read specific entries in a text file
#     y_values_inp.append(float(a[i][38:46]))  # makes sure I backtrack far enough to get all digits(i.e.38 instead of 42)
#     z_values_inp.append(float(a[i][46:54]))
#     positions_inp.append([x_values_inp[i - lines_of_text], y_values_inp[i - lines_of_text], z_values_inp[i - lines_of_text]])
#     identity.append(a[i][12:16])  # hold name of atom (C, H or O)

# Now rotate plane to align with xy plane

# Arrays to hold x,y,z values of each point of interest
# [pt1, pt2, pt3]
plane_x = np.zeros((3, 1))
plane_y = np.zeros((3, 1))
plane_z = np.zeros((3, 1))
# This loop only works because of the way the atoms are spaced in the coordinate file. I am looking at atoms C, C2 and
# C4. Theoretically this will work with any three atoms but I am trying to align the plane of benzene

# if t == 'HII':
#     name = 'HII'
#     for i in range(0, 3):
#         plane_x[i] = float(a[lines_of_text + 2 * i][26:38])
#         plane_y[i] = float(a[lines_of_text + 2 * i][38:46])
#         plane_z[i] = float(a[lines_of_text + 2 * i][46:54])

if t == 'HII':
    name = 'HII'
    plane_atoms = ['C', 'C2', 'C4']
    count = 0
    for i in range(no_atoms):
        if identity[i] in plane_atoms:
            plane_x[count] = float(xyz[0, i])
            plane_y[count] = float(xyz[1, i])
            plane_z[count] = float(xyz[2, i])
            count += 1

if t == 'BCC':
    name = 'MOL'
    plane_x[0] = float(a[lines_of_text + 21][26:38])
    plane_y[0] = float(a[lines_of_text + 21][38:46])
    plane_z[0] = float(a[lines_of_text + 21][46:54])

    plane_x[1] = float(a[lines_of_text + 28][26:38])
    plane_y[1] = float(a[lines_of_text + 28][38:46])
    plane_z[1] = float(a[lines_of_text + 28][46:54])

    plane_x[2] = float(a[lines_of_text + 24][26:38])
    plane_y[2] = float(a[lines_of_text + 24][38:46])
    plane_z[2] = float(a[lines_of_text + 24][46:54])

# vector pointing from point 1 to point 2
v12 = [float(plane_x[1]-plane_x[0]), float(plane_y[1]-plane_y[0]), float(plane_z[1]-plane_z[0])]
v13 = [float(plane_x[2]-plane_x[0]), float(plane_y[2]-plane_y[0]), float(plane_z[2]-plane_z[0])]  # same idea as above

# The cross product of v12 and v13 give a vector that is perpendicular to the plane:
N = np.cross(v12, v13)

N_desired = [0, 0, 1]  # We want the normal vector to be perpendicular to a horizontal plane at the origin

RotationAxis = np.cross(N, N_desired)
theta = math.acos(np.dot(N, N_desired)/(np.linalg.norm(N)*np.linalg.norm(N_desired)))  #  Rotation Angle (radians)

L = [RotationAxis[0]/np.linalg.norm(RotationAxis), RotationAxis[1]/np.linalg.norm(RotationAxis),
                   RotationAxis[2]/np.linalg.norm(RotationAxis)]  # normalized Rotation Axis

# to avoid mistakes: L = [u, v, w]

u = L[0]
v = L[1]
w = L[2]

# Rotation Matrix to rotate a plane:

R = np.zeros((4, 4))
R[3, 3] = 1
R[0, 0] = u**2 + (v**2 + w**2)*math.cos(theta)  # math.cos takes theta in radians by default
R[0, 1] = u*v*(1 - math.cos(theta)) - w*math.sin(theta)
R[0, 2] = u*w*(1 - math.cos(theta)) + v*math.sin(theta)
R[1, 0] = u*v*(1 - math.cos(theta)) + w*math.sin(theta)
R[1, 1] = v**2 + (u**2 + w**2)*math.cos(theta)
R[1, 2] = v*w*(1 - math.cos(theta)) - u*math.sin(theta)
R[2, 0] = u*w*(1 - math.cos(theta)) - v*math.sin(theta)  # math.cos takes theta in radians by default
R[2, 1] = w*v*(1 - math.cos(theta)) + u*math.sin(theta)
R[2, 2] = w**2 + (u**2 + v**2)*math.cos(theta)

# for i in range(0, len(positions_inp)):
#     positions_inp[i].append(1)
#     x = np.dot(R, np.array(positions_inp[i]))
#     positions_inp[i] = [x[0], x[1], x[2]]

b = np.ones([1])
for i in range(np.shape(xyz)[1]):
    coord = np.concatenate((xyz[:, i], b))
    x = np.dot(R, coord)
    xyz[:, i] = x[:3]

# Now translate the structure to the origin
# matrix to translate molecule to origin based on the position of atom 10 (Carbonyl carbon coming off benzene)

# if t == 'HII':
#     translation = np.matrix([[1, 0, 0,-(positions_inp[9][0])], [0, 1, 0,-(positions_inp[9][1])],
#                          [0, 0, 1, -(positions_inp[9][2])], [0, 0, 0, 1]])
if t == 'HII':
    translation = np.matrix([[1, 0, 0,-xyz[0, 9]], [0, 1, 0,-xyz[1, 9]],
                         [0, 0, 1, -xyz[2, 9]], [0, 0, 0, 1]])

# if t == 'BCC':
#     translation = np.matrix([[1, 0, 0, -(positions_inp[25][0])], [0, 1, 0, -(positions_inp[25][1])],
#                          [0, 0, 1, -(positions_inp[25][2])], [0, 0, 0, 1]])
if t == 'BCC':
    translation = np.matrix([[1, 0, 0,-xyz[0, 25]], [0, 1, 0,-xyz[1, 25]],
                         [0, 0, 1, -xyz[2, 25]], [0, 0, 0, 1]])

# for i in range(0, len(positions_inp)):
#     positions_inp[i].append(1)
#     x = np.dot(translation, np.array(positions_inp[i]))
#     positions_inp[i] = [x[0, 0], x[0, 1], x[0, 2]]

b = np.ones([1])
for i in range(np.shape(xyz)[1]):
    coord = np.concatenate((xyz[:, i], b))
    x = np.dot(translation, coord)
    xyz[:, i] = x[0, :3]


# Now rotate the xy coordinates so that the molecule is pointing towards the origin

# if t == 'HII':
#     pt1 = [positions_inp[0][0], positions_inp[0][1]]  # location of C
#     pt2 = [positions_inp[3][0], positions_inp[3][1]]  # location of C3
#     pt3 = [positions_inp[9][0], positions_inp[9][1]]  # location of carbonyl carbon

if t == 'HII':
    pt1 = [xyz[0, 0], xyz[1, 0]]  # location of C
    pt2 = [xyz[0, 3], xyz[1, 3]]  # location of C3
    pt3 = [xyz[0, 9], xyz[1, 9]]  # location of carbonyl carbon

# if t == 'BCC':
#     pt1 = [positions_inp[21][0], positions_inp[21][1]]
#     pt2 = [positions_inp[28][0], positions_inp[28][1]]
#     pt3 = [positions_inp[24][0], positions_inp[24][1]]
if t == 'BCC':
    pt1 = [xyz[0, 21], xyz[1, 21]]
    pt2 = [xyz[0, 28], xyz[1, 28]]
    pt3 = [xyz[0, 24], xyz[1, 24]]

origin = [0, 0]

# find slope between two points

def slope(pt1, pt2):
    m = (pt1[1] - pt2[1])/(pt1[0] - pt2[0])  # slope
    return m

m1 = slope(pt1, pt2)

m2 = 0  # slope of line y = 0

# find angle between lines

theta = -math.atan((m1 - m2)/(1 + m1*m2))

# Find out which quadrant of the xy plane the monomer is sitting in

#   II    |    I
#         |
#   -------------
#         |
#   III   |    IV


def quadrant(pt):  # looks at [x,y] values and determines which quadrant the point is in
    if pt[0] > 0 and pt[1] > 0:
        return 1
    elif pt[0] < 0 and pt[1] < 0:
        return 3
    elif pt[1] < 0 < pt[0]:
        return 4
    elif pt[0] < 0 < pt[1]:
        return 2
    else:
        return 0  # the case where the point lies on the x or y axis

# C and C3 are assumed to be on a line so they should lie in the same quadrant
# Translation matrix
# figure out in which direction the coordinates will be shifted. They are always shift away from the origin

if quadrant(pt1) == 1:  # e.g. in quadrant 1, the x's are shifted in the positive x and positive y directions
    vx = 1
    vy = 1
elif quadrant(pt1) == 2:  # in quadrant 2, the x's are shifted negative and the y's are shifted positive
    vx = -1
    vy = 1
elif quadrant(pt1) == 3:  # in quadrant 3, the x's and y's are shifted down
    vx = -1
    vy = -1
elif quadrant(pt1) == 4:  # in quadrant 4, the x's are shifted positive and the y's are shifted negative
    vx = 1
    vy = -1
# These next three conditionals are very unlikely but are included for completeness and to avoid future errors
elif quadrant(pt1) == 0 and pt1[0] == 0:  # i.e., it lies on the y - axis
    if pt1[1] > 0:  # the point is on the positive y-axis
        vx = 0  # no x-shift
        vy = 1  # shift in the positive y direction
    if pt1[1] < 0:  # the point is on the negative y-axis
        vx = 0  # no x-shift
        vy = -1  # shift in the negative y direction
elif quadrant(pt1) == 0 and pt1[1] == 0:  # i.e., it lies on the x - axis
    if pt1[0] > 0:  # the point is on the positive x-axis
        vx = 1  # shift in the positive x direction
        vy = 0  # no y-shift
    if pt1[0] < 0:  # the point is on the negative y-axis
        vx = -1  # shift in the negative x direction
        vy = 0  # no y-shift

translation = np.matrix([[1, 0, 0, vx*pore_radius*math.cos(theta)], [0, 1, 0, vy*pore_radius*math.sin(theta)],\
                         [0, 0, 1, 0], [0, 0, 0, 1]])

# Apply matrix to all points
# for i in range(0, len(positions_inp)):
#     positions_inp[i].append(1)
#     x = np.dot(translation, np.array(positions_inp[i]))
#     positions_inp[i] = [x[0, 0], x[0, 1], x[0, 2]]
b = np.ones([1])
for i in range(np.shape(xyz)[1]):
    coord = np.concatenate((xyz[:, i], b))
    x = np.dot(translation, coord)
    xyz[:, i] = x[0, :3]

positions = []

# for i in range(0, len(set(layer_distribution))):  # add a list to positions for each unique value of monomers per layer
#     positions.append([])
# for i in range(0, len(positions)):
#     for j in range(0, sorted(list(set(layer_distribution)))[i]):
#         positions[i].append([])

for i in range(len(np.unique(layer_distribution))):  # add a list to positions for each unique value of monomers per layer
    positions.append([])
for i in range(len(positions)):
    for j in range(0, np.sort(np.unique(layer_distribution))[i]):
        positions[i].append([])

x_values = []  # will hold x values in the order that they appear in the positions matrix
y_values = []  # will hold y values in the order that they appear in the positions matrix
z_values = []  # will hold z values in the order that they appear in the positions matrix

# rotate coordinates and store each rotated coordinate as a separate list:

# for j in range(0, len(positions)):
#     for k in range(0, len(positions_inp)):
#         x = np.array(positions_inp[k])
#         for i in range(1, len(positions[j]) + 1):
#             theta = i * math.pi / (len(positions[j]) / 2.0)  # angle to rotate about axis determined from no of monomers per layer
#             # Creates a rotation matrix to rotate input vector about y-axis making a new coordinate at evenly spaced points.
#             # Each rotation belongs to a different monomer's position.The no of points corresponds to the number of monomers
#             Rx = np.zeros((3, 3))  # makes a 3 x 3 zero matrix
#             Rx[0, 0] = math.cos(theta)  # This line and subsequent edits to Rx fills in entries needed for rotation matrix
#             Rx[1, 0] = math.sin(theta)
#             Rx[0, 1] = -math.sin(theta)
#             Rx[1, 1] = math.cos(theta)
#             Rx[2, 2] = 1
#             rot = np.dot(Rx, x)  # multiplies atomic coordinates by the rotation vector to generate new coordinates
#             rot = [float(rot[0]), float(rot[1]), float(rot[2])]  # converts matrix to a list of floats
#             positions[j][i - 1].append(rot)  # appends the atomic coordinates to 'positions'
for j in range(0, len(positions)):
    for k in range(np.shape(xyz)[1]):
        x = np.array(xyz[:, k])
        for i in range(1, len(positions[j]) + 1):
            theta = i * math.pi / (len(positions[j]) / 2.0)  # angle to rotate about axis determined from no of monomers per layer
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
            positions[j][i - 1].append(rot)  # appends the atomic coordinates to 'positions'

f = open('%s' % args.out, 'w')

f.write('This is a .gro file\n')
f.write('%s\n' % sys_atoms)

atom_count = 1
monomer_count = 0
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
        layer_mons = layer_distribution[l*no_layers + k]
        positions_index = sorted(set(layer_distribution)).index(layer_mons)  # The index in positions which should be read from
        for j in range(0, layer_mons):  # iterates over each monomer to create coordinates
            monomer_count += 1
            for i in range(0, no_atoms - no_ions):
                x_values.append(b*dist_bw + positions[positions_index][j][i][0])
                y_values.append(c*dist_bw + positions[positions_index][j][i][1])
                z_values.append(k*dist + positions[positions_index][j][i][2])
                hundreds = int(math.floor(atom_count/100000))
                f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(monomer_count, name, identity[i],
                    atom_count - hundreds*100000, x_values[atom_count - 1]/10.0, y_values[atom_count - 1]/10.0,
                    z_values[atom_count - 1]/10.0) + "\n")
                atom_count += 1

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
        layer_mons = layer_distribution[l*no_layers + k]
        positions_index = sorted(set(layer_distribution)).index(layer_mons)  # The index in positions which should be read from
        for j in range(0, layer_mons):  # iterates over each monomer to create coordinates
            for i in range(0, no_ions):
                x_values.append(b*dist_bw + positions[positions_index][j][no_atoms - (i + 1)][0])
                y_values.append(c*dist_bw + positions[positions_index][j][no_atoms - (i + 1)][1])
                z_values.append(k*dist + positions[positions_index][j][no_atoms - (i + 1)][2])

count = 0
for l in range(0, no_pores):
    for k in range(0, no_layers):
        layer_mons = layer_distribution[l*no_layers + k]
        positions_index = sorted(set(layer_distribution)).index(layer_mons)  # The index in positions which should be read from
        for j in range(0, layer_mons):  # iterates over each monomer to create coordinates
            for i in range(0, no_ions):
                count += 1
                monomer_count += 1  # calling each sodium ion a 'monomer' for ease of numbering
                if atom_count < 100000:
                    hundreds = int(math.floor(atom_count/100000))
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(monomer_count,
                        identity[no_atoms - (no_ions - i)], identity[no_atoms - (no_ions - i)],
                        atom_count - hundreds*100000, x_values[atom_count - 1]/10.0, y_values[atom_count - 1]/10.0,
                        z_values[atom_count - 1]/10.0) + "\n")
                atom_count += 1

f.write('   0.00000   0.00000  0.00000')
f.close()
