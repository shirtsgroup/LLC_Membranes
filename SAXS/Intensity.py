# Calculate intensity as a function of q

import math
import numpy as np
import matplotlib.pyplot as plt
import copy
import random

# Generate qhk for different miller indices in two dimensions

#             -------------
#           /            /
#       b  /            /
#         /            /
#        /            /
#        -------------
#            a

f = open("monomer10_traj.gro", "r")  # .gro file whose positions of Na ions will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

###################     MAKE SURE TO EDIT THE PARAMETERS IN THIS SECTION ACCORDING TO THE SYSTEM     ###################
length_of_simulation = 5000  # picoseconds
no_atoms = 137  # number of atoms in a single monomer
no_layers = 20  # number of layers in the membrane structure
mon_per_layer = 6  # number of monomers in each layer
no_pores = 4  # number of pores
# no_monomers =  # in case there are random distributions of monomers in each layer
no_ion = no_layers*mon_per_layer*no_pores
no_monomers = no_ion  # number of monomers (= number of ions since there is one monomer per ion in the case of sodium)
traj_points = int(len(a) / ((no_atoms + 1)*no_monomers))  # rounds down to the nearest integer. Won't be a problem
# unless there are thousands of trajectory points in which case it'll start rounding up. In that case you can just
# look at the associated .mdp file for the value of traj_points. This is here to be a time saver for the most part
traj_start = 0  # which trajectory to begin analysis at (may need to wait for system to equilibrate)
top_lines = 2  # number of lines at top of each frame
bottom_lines = 2  # number of lines at bottom of each frame

x = []  # list to hold x positions of ions
y = []  # list to hold y positions of ions
z = []  # list to hold z positions of ions. Z axis runs parallel to pore
no_lines_frame = no_atoms*no_layers*mon_per_layer*no_pores + no_ion + bottom_lines  # number of atoms per frame
NA_start = no_atoms*no_layers*mon_per_layer*no_pores + top_lines  # This is the line number where the sodium ions will be listed

# adds all sodium coordinates to respective x, y and z lists
for k in range(0, traj_points):
    for i in range(k*(NA_start + 1 + no_ion) + NA_start, k*(NA_start + 1 + no_ion) + NA_start + no_ion):
        x.append(float(a[i][20:28]))
        y.append(float(a[i][28:36]))
        z.append(float(a[i][36:44]))

# Find the central axis in each pore:
# When looking at a bunch of trajectories, we need to calculate the central axis again each time since it shifts after
# each time step

sum_x_traj = []
sum_y_traj = []
for i in range(0, no_pores):
    sum_x_traj.append([])
    sum_y_traj.append([])

for j in range(traj_start, traj_points):
    for k in range(0, no_pores):  # calculates average x and y values in each pore. Taken to be the center of the pore
        sum_x = 0
        sum_y = 0
        for i in range(j*no_ion + (k * no_ion/4), j*no_ion + ((k + 1) * no_ion/4)):
            sum_x += x[i]
            sum_y += y[i]
        sum_x_traj[k].append(sum_x)
        sum_y_traj[k].append(sum_y)

# average positions to find the location of the central axis in each pore at each frame:
# list structure: [[Pore 1 locations], [Pore 2 locations], [Pore 3 locations], [Pore 4 locations]]

x_axis = []
y_axis = []
for i in range(0, no_pores):
    x_axis.append([])
    y_axis.append([])

constant = 1/float(no_layers*mon_per_layer)  # used to find average. The value of the constant is 1 over the total
# number of molecules in a pore (since there is one of each atom type in each monomer)

for k in range(0, no_pores):  # loops through each pore
    listx = [x_ for x_ in sum_x_traj[k]]  # convert list at the kth index to an array that I can do math on
    listy = [y_ for y_ in sum_y_traj[k]]  # same as line above
    x_axis[k] = np.dot(constant, listx)  # takes average of x_axis location values
    y_axis[k] = np.dot(constant, listy)  # takes average of y_axis location values

# find distance between pores

pore12 = np.zeros(((traj_points - traj_start), 1))
pore13 = np.zeros(((traj_points - traj_start), 1))
pore34 = np.zeros(((traj_points - traj_start), 1))
pore42 = np.zeros(((traj_points - traj_start), 1))
pore14 = np.zeros(((traj_points - traj_start), 1))
pore23 = np.zeros(((traj_points - traj_start), 1))


for i in range(0, traj_points - traj_start):
    pore12[i] = math.sqrt((x_axis[0][i] - x_axis[1][i])**2 + (y_axis[0][i] - y_axis[1][i])**2)
    pore13[i] = math.sqrt((x_axis[0][i] - x_axis[2][i])**2 + (y_axis[0][i] - y_axis[2][i])**2)
    pore34[i] = math.sqrt((x_axis[2][i] - x_axis[3][i])**2 + (y_axis[2][i] - y_axis[3][i])**2)
    pore42[i] = math.sqrt((x_axis[1][i] - x_axis[3][i])**2 + (y_axis[1][i] - y_axis[3][i])**2)
    pore14[i] = math.sqrt((x_axis[0][i] - x_axis[3][i])**2 + (y_axis[0][i] - y_axis[3][i])**2)
    pore23[i] = math.sqrt((x_axis[1][i] - x_axis[2][i])**2 + (y_axis[1][i] - y_axis[2][i])**2)

# Now find average distance between pores

# Pore_list = [pore13, pore14, pore23, pore34, pore12, pore42]
# Sum_list = []
# for i in range(0, len(Pore_list)):
#     Sum_list.append([])
#     Sum_list[i] = 0
#
# for i in range(0, len(Sum_list)):
#     for k in range(traj_start, len(pore12)):
#         Sum_list[i] += Pore_list[i][k]
#
# frames = float(1 / float(traj_points - traj_start))
# Averages = np.dot(Sum_list, frames).tolist()
# Averages.sort()  # Sort values in list from smallest to largest
#
# All of the previous steps just to find this!:
# a = np.mean(Averages[0:5])  # Average of the first 5 points. In a perfect hexagonal lattice all of the values in
# 'Averages' would be the same except for the last entry which is 2*sqrt(3) times every other value in the array
# simply due to the geometry of a hexagon.

# Value of 'a' at each trajectory point
a = []
for i in range(0, traj_points):
    Pore_list = [pore13[i], pore14[i], pore23[i], pore34[i], pore12[i], pore42[i]]
    Pore_list.sort()
    a.append(np.mean(Pore_list[0:5]))


y_axis_reorg = []
x_axis_reorg = []
for i in range(0, traj_points):
    y_axis_reorg.append([y_axis[0][i], y_axis[1][i], y_axis[2][i], y_axis[3][i]])
    x_axis_reorg.append([x_axis[0][i], x_axis[1][i], x_axis[2][i], x_axis[3][i]])

y_unsorted = copy.deepcopy(y_axis_reorg)

for i in range(0, traj_points):
    y_axis_reorg[i].sort()  # sort y_axis values from low to high for each trajectory

# Pair up the correct x values with the sorted y values
x_axis_sorted = []
for k in range(0, traj_points):
    x_axis_sorted.append([])
    for i in range(0, len(y_axis_reorg[k])):
        index = y_axis_reorg[k].index(y_unsorted[k][i])  # find the index of where the y_axis value used to be stored
        x_axis_sorted[k].append(x_axis_reorg[k][index])  # match up the x coordinate corresponding to the y axis coord.

y_axis_sorted = y_axis_reorg  # to help lower confusion

x_shift = np.zeros((no_ion, traj_points))  # will hold positions of shifted x values for each trajectory
y_shift = np.zeros((no_ion, traj_points))  # will hold positions of shifted y values for each trajectory
for i in range(0, traj_points):  # loop through all trajectory point
    for j in range(0, no_pores):  # loop through all pores
        for k in range(0, (no_ion/no_pores)):  # loop through number of ions in each pore
            xshift = min(x_axis[0][i], x_axis[1][i], x_axis[2][i], x_axis[3][i])  # shift by the lowest x_axis value in each trajectory
            yshift = max(y_axis[0][i], y_axis[1][i], y_axis[2][i], y_axis[3][i])  # shift by the highest y_axis value of each trajectory
            x_shift[j*(no_ion/no_pores) + k, i] = x[no_ion*i + j*(no_ion/no_pores) + k] - xshift  # shifting
            y_shift[j*(no_ion/no_pores) + k, i] = y[no_ion*i + j*(no_ion/no_pores) + k] - yshift

x_axis_shifted = []
y_axis_shifted = []
for i in range(0, traj_points):
    x_axis_shifted.append([])
    y_axis_shifted.append([])
    xshift = min(x_axis_sorted[i])
    yshift = max(y_axis_sorted[i])
    for j in range(0, len(x_axis_sorted[i])):
        x_axis_shifted[i].append(x_axis_sorted[i][j] - xshift)
        y_axis_shifted[i].append(y_axis_sorted[i][j] - yshift)

# Rotate All points so that bottom of box is parallel to the x - axis
# Find angle between x axis and bottom line between pores
x_axis_vect = np.array([1, 0])  # unit vector along the x-axis
bottom_vectors = []
angles = []
for i in range(0, traj_points):
    bottom_vectors.append(np.array([(x_axis_shifted[i][2] - x_axis_shifted[i][3]), (y_axis_shifted[i][2] - y_axis_shifted[i][3])]))
    theta = math.acos(np.dot(x_axis_vect, bottom_vectors[i])/(np.linalg.norm(x_axis_vect)*np.linalg.norm(bottom_vectors[i])))
    angles.append(theta)


# Rotation matrix rotates all points counterclockwise by specified angles

x_rot = np.zeros((no_ion, traj_points))  # will hold positions of rotated x values for each trajectory
y_rot = np.zeros((no_ion, traj_points))  # will hold positions of rotated y values for each trajectory

for j in range(0, traj_points):
    R = np.zeros((2, 2))
    R[0, 0] = math.cos(angles[j])
    R[0, 1] = -math.sin(angles[j])
    R[1, 0] = math.sin(angles[j])
    R[1, 1] = math.cos(angles[j])

    for i in range(0, no_ion):
        unrotated = np.matrix([x_shift[i][j], y_shift[i][j]]).T
        rot = np.dot(R, unrotated)
        x_rot[i, j] = rot[0]
        y_rot[i, j] = rot[1]


x_axis_rot = []
y_axis_rot = []
for j in range(0, traj_points):
    x_axis_rot.append([])
    y_axis_rot.append([])
    R = np.zeros((2, 2))
    R[0, 0] = math.cos(angles[j])
    R[0, 1] = -math.sin(angles[j])
    R[1, 0] = math.sin(angles[j])
    R[1, 1] = math.cos(angles[j])
    for i in range(0, no_pores):
        unrotated = np.matrix([x_axis_shifted[j][i], y_axis_shifted[j][i]]).T
        rot = np.dot(R, unrotated)
        x_axis_rot[j].append(rot[0])
        y_axis_rot[j].append(rot[1])

# Convert to fractional coordinates

# First we need the angle characteristic of the unit cell in each frame. They should be about 120, but there are some
# slight fluctuations from frame - to - frame, so it is necessary to calculate this angle for each frame or else the
# fractional coordinates get messed up

frangle = np.zeros((1, traj_points))
for i in range(0, traj_points):
    d10 = np.array((float(x_axis_rot[i][1] - x_axis_rot[i][0]), float(y_axis_rot[i][1] - y_axis_rot[i][0])))
    d13 = np.array((float(x_axis_rot[i][1] - x_axis_rot[i][3]), float(y_axis_rot[i][1] - y_axis_rot[i][3])))
    dot = np.dot(d10, d13)
    frangle[0, i] = math.acos(dot/(np.linalg.norm(d10)*np.linalg.norm(d13)))

x_frac = np.zeros((no_ion, traj_points))  # will hold fractional coordinates of x values for each trajectory frame
y_frac = np.zeros((no_ion, traj_points))  # will hold fractional coordinates of y values for each trajectory frame
for k in range(0, traj_points):
    for i in range(0, len(x_rot)):
        x_frac[i, k] = (1/a[k])*x_rot[i, k] - y_rot[i, k]*math.cos(frangle[0, k])/(a[k]*math.sin(frangle[0, k]))
        y_frac[i, k] = y_rot[i, k]/(a[k]*math.sin(frangle[0, k]))

# Now take care of the points that lie outside of the box enclosed by [0,0], [1, 0], [0, 1], [1, 1] (enforce PBC's)
x_frac_box = np.zeros((no_ion, traj_points))  # will hold fractional coordinates of x values for each trajectory frame
y_frac_box = np.zeros((no_ion, traj_points))  # will hold fractional coordinates of y values for each trajectory frame
for k in range(0, traj_points):
    for i in range(0, no_ion):
        if x_frac[i, k] < 0 and y_frac[i, k] > 0:
            x_frac_box[i, k] = x_frac[i, k] + 1
            y_frac_box[i, k] = y_frac[i, k] - 1
        elif 1 >= x_frac[i, k] >= 0 and y_frac[i, k] > 0:
            x_frac_box[i, k] = x_frac[i, k]
            y_frac_box[i, k] = y_frac[i, k] - 1
        elif x_frac[i, k] < 0 and -1 <= y_frac[i, k] <= 0:
            x_frac_box[i, k] = x_frac[i, k] + 1
            y_frac_box[i, k] = y_frac[i, k]
        elif x_frac[i, k] < 0 and y_frac[i, k] < - 1:
            x_frac_box[i, k] = x_frac[i, k] + 1
            y_frac_box[i, k] = y_frac[i, k] + 1
        elif 1 >= x_frac[i, k] >= 0 and y_frac[i, k] < -1:
            x_frac_box[i, k] = x_frac[i, k]
            y_frac_box[i, k] = y_frac[i, k] + 1
        elif x_frac[i, k] > 1 and y_frac[i, k] < -1:
            x_frac_box[i, k] = x_frac[i, k] - 1
            y_frac_box[i, k] = y_frac[i, k] + 1
        elif x_frac[i, k] > 1 and y_frac[i, k] > 0:
            x_frac_box[i, k] = x_frac[i, k] - 1
            y_frac_box[i, k] = y_frac[i, k] - 1
        elif x_frac[i, k] > 1 and -1 <= y_frac[i, k] <= 0:
            x_frac_box[i, k] = x_frac[i, k] - 1
            y_frac_box[i, k] = y_frac[i, k]
        else:
            y_frac_box[i, k] = y_frac[i, k]
            x_frac_box[i, k] = x_frac[i, k]

plt.plot(x_frac_box, y_frac_box)
plt.show()
# re = 2.81794*10**-15  # radius of the electron [m]
# Sc = a**2  # area of 2D unit cell [nm^2]
# L = 28000  # membrane thickness [nm]

# def qhk(h, k, a):
#     # h, k are 2D miller indices
#     # a is the box vector (assumed equal to b)
#     return 4*math.pi/(math.sqrt(3)*a)*math.sqrt(h**2 + k**2 + h*k)  # based on supplementary info from imperor-clerc


# def scattering_factor(a, b, c, d):
#     # d: 1/2d where d is the d-spacing
#     f = 0  # intialize summation
#     for i in range(0, 4):  # for a1, b1, a2, b2, a3, b3, a4, b4
#         f += a[i]*math.exp(-b[i]*d**2)  # this function was used to fit the data in the Crystallographic tables
#     f += c  # the constant c is contained outside the summation
#     return f
#
#
# def structure_factor(h, k, a, x, y):   # |Fhk|^2
#     d = a/math.sqrt(h**2 + k**2)  # d based on miller indices of this plane defined by h and k
#     ratio = 1/(2*d)  # sin(theta) / lambda = 1/ 2d
#     f = scattering_factor(a_p, b_p, c_p, ratio)  # find scattering factor for this value of 1/2d
#     sum1 = 0
#     sum2 = 0
#     for i in range(0, len(x)):
#         sum1 += f*math.cos(2*math.pi*(h*x[i] + k*y[i]))
#         sum2 += f*math.sin(2*math.pi*(h*x[i] + k*y[i]))
#     return sum1**2 + sum2**2
#
# x_rand = []
# y_rand = []
# n_rand = 480
# for i in range(0, n_rand):
#     x_rand.append(random.random())
#     y_rand.append(random.random())
#
# Intensities = []
# Trajectory = np.linspace(0, traj_points, traj_points)
# for i in range(0, len(Trajectory)):
#     Intensities.append(structure_factor(2, 0, a[i], x_frac_box[:, i], y_frac_box[:, i]))
#     # Intensities.append(structure_factor(1, 0, 4, x_rand, y_rand))
#
# print np.mean(Intensities[30:49])
# plt.plot(Trajectory, Intensities)
# plt.show()
# def Laue(x, N):  # Laue Function
#     # x is the q vector
#     # N is qhk
#     return (math.sin(N*x) / math.sin(x))**2
#
#
# def Intensity(N):  # evaluate I as a function of q for different values of N
#     q_values = np.zeros((N, N))  # matrix to hold values of different qhk
#
#     # multiplicities of each miller plane:
#     m = 12*np.ones((N, N))  # matrix holding multiplicity of each miller index combination
#     # if h does not equal k and k is not zero, the multiplicity is 12
#     m[:, 0] = 6  # for all h0, multiplicity = 6
#     for i in range(0, N - 1):  # for all h = k, multiplicity = 6
#         m[i, i + 1] = 6
#
#     # Calculate all qhk as defined by the size of N
#     for i in range(0, N):
#         for k in range(0, N):
#         #for k in range(0, i + 1): # only calculates main diagonal and values below since q12 = q21. Saves comp time
#             q_values[i, k] = qhk(i + 1, k, a)
#
#     q = np.linspace(0.1, 0.5, 100)  # plot enough points to make a continuous looking function
#     I_L = np.zeros((len(q), 1))
#
#     for i in range(0, len(q)):  # evaluate all of the different q
#         Laue_eval = np.zeros((N, N))  # intialize matrix to hold evaluations of Laue function
#         Summation = np.zeros((N, N))  # functional evaluations of the Laue function multiplied by their multiplicity
#         Fhk = np.zeros((N, N))  # evaluate |Fhk|^2 term for each miller index we are looking at
#         for j in range(0, N):  # go through all rows
#             for k in range(0, j + 1):  # go through all columns (not looking at values above main diagonal)
#                 Fhk[j, k] = structure_factor(j + 1, k, a)
#                 Laue_eval[j, k] = Laue(q[i], q_values[j, k])  # evaluate laue function using various q and qhk values
#                 Summation[j, k] = Laue_eval[j, k] * m[j, k] * Fhk[j, k]  # multiply by its multiplicity
#         I_L[i] = (re**2*2*math.pi**2/q[i]**2)*sum(sum(Summation))
#
#     return I_L, q
#
#
# I_L, q = Intensity(1)
# plt.plot(q, I_L)
# plt.title('Sum of Laue functional evaluations')
# plt.xlabel('q [nm^-1]')
# plt.ylabel('Intensity')
# plt.show()
