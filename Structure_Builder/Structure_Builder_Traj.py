import os
import numpy as np
import math
from random import randint
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
dist = 10  # distance between layers (units tbd)
lines_of_text = 4


files = []  # list to hold the names of files
for i in os.listdir('/home/bcoscia/PycharmProjects/LLC_Structure_Builder/Monomer_Configurations'):  # looks in folder
    # located in this path
    if i.endswith(".pdb"):  # if the file is a .gro it will be taken.
        files.append(i)  # adds file to list f


a = []
for i in range(0, len(files)):
    b = []
    f = open('Monomer_Configurations/%s' %files[i], "r")
    for line in f:
        b.append(line)
    a.append(b)

# Make a list to hold the identity of atoms at each position. Only done once since all files are of the same format
identity = []
for i in range(4, 142):
    identity.append(a[0][i][13:16])  # hold name of atom (C, H, NA or O)

x_values_inp = []  # list to hold input values of x stored from .pdb file
y_values_inp = []  # list to hold input values of y stored from .pdb file
z_values_inp = []  # list to hold input values of z stored from .pdb file
positions_inp = []  # holds x, y, z coordinates of input .pdb file
for k in range(0, len(files)):
    x_values_file = []  # list to hold input values of x stored from .pdb file for this iteration
    y_values_file = []  # list to hold input values of y stored from .pdb file for this iteration
    z_values_file = []  # list to hold input values of z stored from .pdb file for this iteration
    positions_file = []  # holds x, y, z coordinates of input .pdb file for this iteration
    # read specific entries in a text file
    for i in range(4, 142):  # searches relevant lines of text in file, f, being read
        # There are 137 atoms excluding sodium
        x_values_file.append(float(a[k][i][26:38]))  # Use this to read specific entries in a text file
        y_values_file.append(float(a[k][i][38:46]))  # makes sure I backtrack far enough to get all digits(i.e.38 instead of 42)
        z_values_file.append(float(a[k][i][46:54]))
        positions_file.append([x_values_file[i - 4], y_values_file[i - 4], z_values_file[i - 4]])  # i - 4 since i starts at 4
    x_values_inp.append(x_values_file)
    y_values_inp.append(y_values_file)
    z_values_inp.append(z_values_file)
    positions_inp.append(positions_file)

# Now rotate each monomer's benzene plane to align with xy plane

# Arrays to hold x,y,z values of each point of interest for each monomer individually
# [pt1, pt2, pt3]
x_benz = []
y_benz = []
z_benz = []
for i in range(0, len(files)):
    x_benz.append([])
    y_benz.append([])
    z_benz.append([])

for j in range(0, len(files)):
    plane_x = np.zeros((3, 1))
    plane_y = np.zeros((3, 1))
    plane_z = np.zeros((3, 1))
# This loop only works because of the way the atoms are spaced in the coordinate file. I am looking at atoms C, C2 and
# C4. Theoretically this will work with any three atoms but I am trying to align the plane of benzene
    for i in range(0, 3):
        plane_x[i] = float(a[j][lines_of_text + 2*i][26:38])
        plane_y[i] = float(a[j][lines_of_text + 2*i][38:46])
        plane_z[i] = float(a[j][lines_of_text + 2*i][46:54])
    x_benz[j] = plane_x
    y_benz[j] = plane_y
    z_benz[j] = plane_z


# vector pointing from point 1 to point 2
v12 = []
v13 = []
for i in range(0, len(files)):
    v12.append([])
    v13.append([])

for i in range(0, len(files)):
    v12[i] = [float(x_benz[i][1]-x_benz[i][0]), float(y_benz[i][1]-y_benz[i][0]), float(z_benz[i][1]-z_benz[i][0])]
    v13[i] = [float(x_benz[i][2]-x_benz[i][0]), float(y_benz[i][2]-y_benz[i][0]), float(z_benz[i][2]-z_benz[i][0])]

# The cross product of v12 and v13 give a vector that is perpendicular to the plane:
N = []
for i in range(0, len(files)):
    N.append([])

for i in range(0, len(files)):
    N[i] = np.cross(v12[i], v13[i])

N_desired = [0, 0, 1]  # We want the normal vector to be perpendicular to a horizontal plane at the origin

RotationAxis = []
theta = []
L = []
for i in range(0, len(files)):
    RotationAxis.append([])
    theta.append([])
    L.append([])

for i in range(0, len(files)):
    RotationAxis[i] = np.cross(N[i], N_desired)
    theta[i] = math.acos(np.dot(N[i], N_desired)/(np.linalg.norm(N[i])*np.linalg.norm(N_desired)))  #  Rotation Angle (radians)
    # normalize the rotation axis vectors
    L[i] = [RotationAxis[i][0]/np.linalg.norm(RotationAxis[i]), RotationAxis[i][1]/np.linalg.norm(RotationAxis[i]),
            RotationAxis[i][2]/np.linalg.norm(RotationAxis[i])]

# to avoid mistakes: L = [u, v, w]
u = []
v = []
w = []
for i in range(0, len(files)):
    u.append(L[i][0])
    v.append(L[i][1])
    w.append(L[i][2])

# Rotation Matrix to rotate a plane:
for j in range(0, len(files)):
    R = np.zeros((4, 4))
    R[3, 3] = 1
    R[0, 0] = u[j]**2 + (v[j]**2 + w[j]**2)*math.cos(theta[j])  # math.cos takes theta in radians by default
    R[0, 1] = u[j]*v[j]*(1 - math.cos(theta[j])) - w[j]*math.sin(theta[j])
    R[0, 2] = u[j]*w[j]*(1 - math.cos(theta[j])) + v[j]*math.sin(theta[j])
    R[1, 0] = u[j]*v[j]*(1 - math.cos(theta[j])) + w[j]*math.sin(theta[j])
    R[1, 1] = v[j]**2 + (u[j]**2 + w[j]**2)*math.cos(theta[j])
    R[1, 2] = v[j]*w[j]*(1 - math.cos(theta[j])) - u[j]*math.sin(theta[j])
    R[2, 0] = u[j]*w[j]*(1 - math.cos(theta[j])) - v[j]*math.sin(theta[j])  # math.cos takes theta in radians by default
    R[2, 1] = w[j]*v[j]*(1 - math.cos(theta[j])) + u[j]*math.sin(theta[j])
    R[2, 2] = w[j]**2 + (u[j]**2 + v[j]**2)*math.cos(theta[j])

    for i in range(0, len(positions_inp[j])):
        positions_inp[j][i].append(1)
        x = np.dot(R, np.array(positions_inp[j][i]))
        positions_inp[j][i] = [x[0], x[1], x[2]]

# Now translate the structure to the origin
# matrix to translate molecule to origin based on the position of atom 10 (Carbonyl carbon coming off benzene)
for i in range(0, len(files)):
    translation = np.matrix([[1, 0, 0, -(positions_inp[i][9][0])], [0, 1, 0, -(positions_inp[i][9][1])],
                             [0, 0, 1, -(positions_inp[i][9][2])], [0, 0, 0, 1]])
    for j in range(0, len(positions_inp)):
        positions_inp[i][j].append(1)
        x = np.dot(translation, np.array(positions_inp[i][j]))
        positions_inp[i][j] = [x[0, 0], x[0, 1], x[0, 2]]

# Now rotate the xy coordinates so that the molecule is pointing towards the origin

pt1 = []
pt2 = []

for i in range(0, len(files)):
    pt1.append([])
    pt2.append([])

for i in range(0, len(files)):
    pt1[i] = [positions_inp[i][0][0], positions_inp[i][0][1]]  # location of C
    pt2[i] = [positions_inp[i][3][0], positions_inp[i][3][1]]  # location of C3

# find slope between two points
def slope(pt1, pt2):
    m = (pt1[1] - pt2[1])/(pt1[0] - pt2[0])  # slope
    return m

m1 = []
m2 = []
angle = []
for i in range(0, len(files)):
    m1.append([])
    m2.append([])
    angle.append([])

for i in range(0, len(files)):
    m1[i] = slope(pt1[i], pt2[i])
    m2[i] = 0  # slope of line y = 0

# find angle between lines
for i in range(0, len(files)):
    angle[i] = -math.atan((m1[i] - m2[i])/(1 + m1[i]*m2[i]))

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
    elif pt[1] < 0 < pt[0] :
        return 4
    elif pt[0] < 0 < pt[1] :
        return 2
    else:
        return 0  # the case where the point lies on the x or y axis

# C and C3 are assumed to be on a line so they should lie in the same quadrant

#  figure out in which direction the coordinates will be shifted. They are always shift away from the origin
for i in range(0, len(files)):
    vx = 0  # case where the monomer is already at the origin
    vy = 0
    if quadrant(pt1[i]) == 1:  # e.g. in quadrant 1, the x's are shifted in the positive x and positive y directions
        vx = 1
        vy = 1
    elif quadrant(pt1[i]) == 2:  # in quadrant 2, the x's are shifted negative and the y's are shifted positive
        vx = -1
        vy = 1
    elif quadrant(pt1[i]) == 3:  # in quadrant 3, the x's and y's are shifted down
        vx = -1
        vy = -1
    elif quadrant(pt1[i]) == 4:  # in quadrant 4, the x's are shifted positive and the y's are shifted negative
        vx = 1
        vy = -1
    # These next three conditionals are very unlikely but are included for completeness and to avoid future errors
    elif quadrant(pt1[i]) == 0 and pt1[i][0] == 0:  # i.e., it lies on the y - axis
        if pt1[i][1] > 0:  # the point is on the positive y-axis
            vx = 0  # no x-shift
            vy = 1  # shift in the positive y direction
        if pt1[i][1] < 0:  # the point is on the negative y-axis
            vx = 0  # no x-shift
            vy = -1  # shift in the negative y direction
    elif quadrant(pt1[i]) == 0 and pt1[i][1] == 0:  # i.e., it lies on the x - axis
        if pt1[i][0] > 0:  # the point is on the positive x-axis
            vx = 1  # shift in the positive x direction
            vy = 0  # no y-shift
        if pt1[i][0] < 0:  # the point is on the negative y-axis
            vx = -1  # shift in the negative x direction
            vy = 0  # no y-shift

    translation = np.matrix([[1, 0, 0, vx*pore_radius*math.cos(angle[i])], [0, 1, 0, vy*pore_radius*math.sin(angle[i])],\
                         [0, 0, 1, 0], [0, 0, 0, 1]])

    # Apply matrix to all points
    for j in range(0, len(positions_inp[i])):
        positions_inp[i][j].append(1)
        x = np.dot(translation, np.array(positions_inp[i][j]))
        positions_inp[i][j] = [x[0, 0], x[0, 1], x[0, 2]]


positions = [[], [], [], [], [], [], [], [], [], []]
x_values = []  # will hold x_values in the order that they appear in the positions matrix
y_values = []  # will hold y_values in the order that they appear in the positions matrix
z_values = []  # will hold z_values in the order that they appear in the positions matrix

# rotate coordinates and store each rotated coordinate as a separate list:

for l in range(0, len(files)):  # will rotate coordinates in all files
    positions_file = [[], [], [], [], [], [], []]
    for k in range(0, len(positions_inp[l])):  # iterates through coordinates of all atoms in each file
        x = np.array(positions_inp[l][k])  # looking at a specific atom's xyz coordinates
        for i in range(1, no_monomers + 1):  # rotates coordinates at equal angles around axis
            theta = i * math.pi / (no_monomers / 2.0)  # angle to rotate about axis determined from no of monomers per layer
            # Creates a rotation matrix to rotate input vector about y-axis making a new coordinate at evenly spaced points.
            # Each rotation belongs to a different monomer's position.The no of points corresponds to the number of monomers
            Rx = np.zeros((3, 3))  # makes a 3 x 3 zero matrix
            Rx[0, 0] = math.cos(theta)  # This line and subsequent edits to Rx creates rotation matrix
            Rx[1, 0] = math.sin(theta)
            Rx[0, 1] = -math.sin(theta)
            Rx[1, 1] = math.cos(theta)
            Rx[2, 2] = 1
            rot = np.dot(Rx, x)  # multiplies atomic coordinates by the rotation vector to generate new coordinates
            rot = [float(rot[0]), float(rot[1]), float(rot[2])]  # converts matrix to a list of floats
            positions_file[i - 1].append(rot)  # appends the atomic coordinates to 'positions'
    positions[l].append(positions_file)

#Now randomly choose monomers from 'positions' to create the final positions of all atoms


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
    for k in range(0, no_layers):  # creates a new layer
        final_positions = [[], [], [], [], [], [], []]
        for g in range(0, no_monomers):
            s = randint(0, 9)  # random number corresponding to which monomer file will be used
            for h in range(0, no_atoms):
                final_positions[g].append(positions[s][0][g][h])
        for j in range(0, no_monomers):  # iterates over each monomer to create coordinates
            for i in range(0, no_atoms - 1):
                x_values.append(b*dist_bw + final_positions[j][i][0])
                y_values.append(c*dist_bw + final_positions[j][i][1])
                z_values.append(k*dist + final_positions[j][i][2])
                print '{:6s}{:5d}  {:4s} {:8d}   {:8.3f} {:8.3f}{:8.3f}  {:4.2f}  {:4.2f}          {:>2}{:2s}'\
                    .format('ATOM', i + 1 + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l,
                    identity[i], 0, x_values[i+(no_atoms - 1)*j+(no_atoms - 1)*no_monomers*k+(no_atoms - 1)*no_monomers*no_layers*l],
                    y_values[i + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l],
                    z_values[i + (no_atoms - 1)*j + (no_atoms - 1)*no_monomers*k + (no_atoms - 1)*no_monomers*no_layers*l], 0.00,
                    0.00, 'C', '+0')

# Now add sodium:


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
            x_values.append(b*dist_bw + final_positions[j][no_atoms - 1][0])
            y_values.append(c*dist_bw + final_positions[j][no_atoms - 1][1])
            z_values.append(k*dist + final_positions[j][no_atoms -1][2])
            print '{:6s}{:5d}  {:4s} {:8d}   {:8.3f} {:8.3f}{:8.3f}  {:4.2f}  {:4.2f}          {:>2}{:2s}'\
                .format('ATOM', 1 + (no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers,
                identity[no_atoms -1], 0, x_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers],
                y_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers],
                z_values[(no_atoms - 1)*no_layers*no_pores*no_monomers + j + k*no_monomers + l*no_layers*no_monomers], 0.00,
                0.00, 'C', '+0')

# c1 = []  # connectivity column 1 from .pdb file
# c2 = []  # connectivity column 2 from .pdb file
# c3 = []  # connectivity column 3 from .pdb file
# c4 = []  # connectivity column 4 from .pdb file
# c5 = []  # connectivity column 5 from .pdb file
#
# for l in range(0, no_pores):
#     for k in range(0, no_layers):
#         for j in range(0, no_monomers):
#             for i in range(143, 203):  # loops through rows which contain connectivity information
#                 if int(a[i][6:11]) != 0:  # need to retain 0's since they mean that there is no connection at that point
#                     c1.append(j*(no_atoms - 1) + int(a[i][6:11]))  # searches for up to a 5 digit number
#                 else:
#                     c1.append(0)
#                 if int(a[i][11:16]) != 0:
#                     c2.append(j*(no_atoms - 1) + int(a[i][11:16]))
#                 else:
#                     c2.append(0)
#                 if int(a[i][16:21]) != 0:
#                     c3.append(j*(no_atoms - 1) + int(a[i][16:21]))
#                 else:
#                     c3.append(0)
#                 if int(a[i][21:26]) != 0:
#                     c4.append(j*(no_atoms - 1) + int(a[i][21:26]))
#                 else:
#                     c4.append(0)
#                 if int(a[i][26:31]) != 0:
#                     c5.append(j*(no_atoms - 1) + int(a[i][26:31]))
#                 else:
#                     c5.append(0)
#                 print '{:6s}{:5d}{:5d}{:5d}{:5d}{:5d}                                         {:4s}{:4d}'\
#                     .format('CONECT', c1[j * 60 + i - 143], c2[j * 60 + i - 143], c3[j * 60 + i - 143],
#                     c4[j * 60 + i - 143], c5[j * 60 + i - 143], 'NONE', j*60 + 5 + (no_atoms -1)*(no_monomers-1) + i + 1)
# # line numbering (last column) is a bit off but it doesn't matter
# print 'END'