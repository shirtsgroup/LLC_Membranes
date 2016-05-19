#!/usr/bin/python

import numpy as np
import math
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import griddata
from matplotlib import animation
import argparse

# Pore Numbering : Each corner is the location of the center of the pore
#      Pore 1 ------------- Pore 2'
#           /            /
#          /            /
#         /            /
#        /            /
# Pore 4 ------------- Pore 3'

parser = argparse.ArgumentParser(description = 'Run Cylindricity script')
parser.add_argument('-i', '--input', default='Monomer1_Traj.pdb', help = 'Path to input file')
args = parser.parse_args()

f = open(args.input, "r")  # .gro file whose positions of Na ions will be read
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
carb_index = 9 + top_lines  # this is where the first carbonyl carbon is listed in the file

# adds all sodium coordinates to respective x, y and z lists
for k in range(0, traj_points):
    for i in range(k*(NA_start + 1 + no_ion) + NA_start, k*(NA_start + 1 + no_ion) + NA_start + no_ion):
        x.append(float(a[i][20:28]))
        y.append(float(a[i][28:36]))
        z.append(float(a[i][36:44]))

x_carb = []  # list to hold x positions of carbonyl carbons
y_carb = []  # list to hold y positions of carbonyl carbons
z_carb = []  # list to hold z positions of carbonyl carbons

# extract coordinates of all carbonyl carbons (C6)
for k in range(0, traj_points):
    for i in range(0, no_monomers):
        index = i*no_atoms + k*(no_lines_frame + 1) + carb_index  # location of C6 as a function of i and k
        x_carb.append(float(a[index][20:28]))
        y_carb.append(float(a[index][28:36]))
        z_carb.append(float(a[index][36:44]))

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


# Do the same with carbonyl carbons

sum_x_traj_carb = []
sum_y_traj_carb = []
for i in range(0, no_pores):
    sum_x_traj_carb.append([])
    sum_y_traj_carb.append([])

for j in range(traj_start, traj_points):
    for k in range(0, no_pores):  # calculates average x and y values in each pore. Taken to be the center of the pore
        sum_x = 0
        sum_y = 0
        for i in range(j*no_ion + (k * no_ion/4), j*no_ion + ((k + 1) * no_ion/4)):
            sum_x += x_carb[i]
            sum_y += y_carb[i]
        sum_x_traj_carb[k].append(sum_x)
        sum_y_traj_carb[k].append(sum_y)

# average positions to find the location of the central axis in each pore at each frame:
# list structure: [[Pore 1 locations], [Pore 2 locations], [Pore 3 locations], [Pore 4 locations]]

x_axis = []
y_axis = []
x_axis_carb = []
y_axis_carb = []
for i in range(0, no_pores):
    x_axis.append([])
    y_axis.append([])
    x_axis_carb.append([])
    y_axis_carb.append([])

constant = 1/float(no_layers*mon_per_layer)  # used to find average. The value of the constant is 1 over the total
# number of molecules in a pore (since there is one of each atom type in each monomer)

for k in range(0, no_pores):  # loops through each pore
    listx = [x_ for x_ in sum_x_traj[k]]  # convert list at the kth index to an array that I can do math on
    listy = [y_ for y_ in sum_y_traj[k]]  # same as line above
    listx_carb = [x_ for x_ in sum_x_traj_carb[k]]  # convert list at the kth index to an array that I can do math on
    listy_carb = [y_ for y_ in sum_y_traj_carb[k]]  # same as line above
    x_axis[k] = np.dot(constant, listx)  # takes average of x_axis location values
    y_axis[k] = np.dot(constant, listy)  # takes average of y_axis location values
    x_axis_carb[k] = np.dot(constant, listx_carb)
    y_axis_carb[k] = np.dot(constant, listy_carb)

# find distance between pores

pore12 = np.zeros(((traj_points - traj_start), 1))
pore13 = np.zeros(((traj_points - traj_start), 1))
pore34 = np.zeros(((traj_points - traj_start), 1))
pore42 = np.zeros(((traj_points - traj_start), 1))
pore14 = np.zeros(((traj_points - traj_start), 1))
pore23 = np.zeros(((traj_points - traj_start), 1))
pore12_carb = np.zeros(((traj_points - traj_start), 1))
pore13_carb = np.zeros(((traj_points - traj_start), 1))
pore34_carb = np.zeros(((traj_points - traj_start), 1))
pore42_carb = np.zeros(((traj_points - traj_start), 1))
pore14_carb = np.zeros(((traj_points - traj_start), 1))
pore23_carb = np.zeros(((traj_points - traj_start), 1))

for i in range(0, traj_points - traj_start):
    pore12[i] = math.sqrt((x_axis[0][i] - x_axis[1][i])**2 + (y_axis[0][i] - y_axis[1][i])**2)
    pore13[i] = math.sqrt((x_axis[0][i] - x_axis[2][i])**2 + (y_axis[0][i] - y_axis[2][i])**2)
    pore34[i] = math.sqrt((x_axis[2][i] - x_axis[3][i])**2 + (y_axis[2][i] - y_axis[3][i])**2)
    pore42[i] = math.sqrt((x_axis[1][i] - x_axis[3][i])**2 + (y_axis[1][i] - y_axis[3][i])**2)
    pore14[i] = math.sqrt((x_axis[0][i] - x_axis[3][i])**2 + (y_axis[0][i] - y_axis[3][i])**2)
    pore23[i] = math.sqrt((x_axis[1][i] - x_axis[2][i])**2 + (y_axis[1][i] - y_axis[2][i])**2)
    pore12_carb[i] = math.sqrt((x_axis_carb[0][i] - x_axis_carb[1][i])**2 + (y_axis_carb[0][i] - y_axis_carb[1][i])**2)
    pore13_carb[i] = math.sqrt((x_axis_carb[0][i] - x_axis_carb[2][i])**2 + (y_axis_carb[0][i] - y_axis_carb[2][i])**2)
    pore34_carb[i] = math.sqrt((x_axis_carb[2][i] - x_axis_carb[3][i])**2 + (y_axis_carb[2][i] - y_axis_carb[3][i])**2)
    pore42_carb[i] = math.sqrt((x_axis_carb[1][i] - x_axis_carb[3][i])**2 + (y_axis_carb[1][i] - y_axis_carb[3][i])**2)
    pore14_carb[i] = math.sqrt((x_axis_carb[0][i] - x_axis_carb[3][i])**2 + (y_axis_carb[0][i] - y_axis_carb[3][i])**2)
    pore23_carb[i] = math.sqrt((x_axis_carb[1][i] - x_axis_carb[2][i])**2 + (y_axis_carb[1][i] - y_axis_carb[2][i])**2)

plt.figure(1)
time_pts = range(0, (traj_points - traj_start))
plt.plot(time_pts, pore12, label='1-2')
plt.plot(time_pts, pore13, label='1-3')
plt.plot(time_pts, pore34, label='3-4')
plt.plot(time_pts, pore42, label='4-2')
plt.plot(time_pts, pore14, label='4-1')
plt.plot(time_pts, pore23, label='2-3')
plt.title('Pore-To-Pore Distance Equilibration')
plt.ylabel('Pore-To-Pore Distance [nm]')
plt.xlabel('Frame')
plt.legend(loc=1)

plt.figure(2)
time_pts = range(0, (traj_points - traj_start))
plt.plot(time_pts, pore12_carb, label='1-2')
plt.plot(time_pts, pore13_carb, label='1-3')
plt.plot(time_pts, pore34_carb, label='3-4')
plt.plot(time_pts, pore42_carb, label='4-2')
plt.plot(time_pts, pore14_carb, label='4-1')
plt.plot(time_pts, pore23_carb, label='2-3')
plt.title('Pore-To-Pore Distance Equilibration based on Carbonyl Carbon')
plt.ylabel('Pore-To-Pore Distance [nm]')
plt.xlabel('Frame')
plt.legend(loc=1)

# Find the distance from central axis
# need to save each value of distance in order to calculate deviation. 4 lists in a list

dist_list = []
dist_list_carb = []
for i in range(0, no_pores):
    dist_list.append([])
    dist_list_carb.append([])

for l in range(0, traj_points - traj_start):
    for k in range(0, no_pores):
        for j in range(0, no_layers):
            for i in range(0, mon_per_layer):  # loops through each pore individually
                dist = math.sqrt((x[i + mon_per_layer*j + mon_per_layer*no_layers*k + l*mon_per_layer*no_layers*no_pores] - x_axis[k][l])**2 +
                                  (y[i + mon_per_layer*j + mon_per_layer*no_layers*k + l*mon_per_layer*no_layers*no_pores] - y_axis[k][l])**2)
                dist_carb = math.sqrt((x_carb[i + mon_per_layer*j + mon_per_layer*no_layers*k + l*mon_per_layer*no_layers*no_pores] - x_axis_carb[k][l])**2 +
                  (y_carb[i + mon_per_layer*j + mon_per_layer*no_layers*k + l*mon_per_layer*no_layers*no_pores] - y_axis_carb[k][l])**2)
                # magnitude of vector pointing from center to point where ion is stationed
                dist_list[k].append(dist)
                dist_list_carb[k].append(dist_carb)

# find averages at each trajectory point to monitor equilibration:

traj_averages = np.zeros(((traj_points - traj_start), 1))  # array to hold average values
traj_averages_carb = np.zeros(((traj_points - traj_start), 1))  # array to hold average values
pts = range(0, (traj_points - traj_start))  # generate points for plotting traj_averages versus trajectory point

for j in range(0, traj_points - traj_start):  # will generate a point for each trajectory point
    traj_sum = 0
    traj_sum_carb = 0
    for i in range(0, no_layers*mon_per_layer):  # loops through all atoms in a pore
        for k in range(0, no_pores):  # there is a different list for each pore in dist_list
            traj_sum += dist_list[k][j*mon_per_layer*no_layers + i]
            traj_sum_carb += dist_list_carb[k][j*mon_per_layer*no_layers + i]
    traj_pt_avg = traj_sum / (no_pores*no_layers*mon_per_layer)  # finds average distance for this trajectory point
    traj_pt_avg_carb = traj_sum_carb / (no_pores*no_layers*mon_per_layer)
    traj_averages[j] = traj_pt_avg
    traj_averages_carb[j] = traj_pt_avg_carb

plt.figure(3)
plt.plot(pts, traj_averages, label = 'Sodium')
plt.plot(pts, traj_averages_carb, label = 'Carbonyl')
plt.title('%s Picosecond Simulation' %length_of_simulation)
plt.ylabel('Average Distance from Pore Center (nm)')
plt.xlabel('Trajectory Point')

# Calculate overall average values and deviation

avg = np.zeros((no_pores, 1))
for i in range(0, no_pores):
    avg[i] = np.mean(dist_list[i])

Average_Distance = np.mean(avg)

std = np.zeros((no_pores, 1))
for i in range(0, no_pores):
    std[i] = np.std(dist_list[i])

Average_Deviation = np.mean(std)

# combine the four lists from dist_list into one list:
dist_one_list = np.zeros((no_pores*(traj_points-traj_start)*no_layers*mon_per_layer))
dist_one_list_carb = np.zeros((no_pores*(traj_points-traj_start)*no_layers*mon_per_layer))
for i in range(0, no_pores):
    for k in range(0, (traj_points-traj_start)*no_layers*mon_per_layer):
        dist_one_list[i*(traj_points-traj_start)*no_layers*mon_per_layer + k] = dist_list[i][k]
        dist_one_list_carb[i*(traj_points-traj_start)*no_layers*mon_per_layer + k] = dist_list_carb[i][k]


# Create a distribution of distances
bar_width = 0.05
distribution_pts = np.arange(0, max(dist_one_list) + bar_width, bar_width)  # list to group distances based on being between two points
distribution_pts_carb = np.arange(0, max(dist_one_list_carb) + bar_width, bar_width)
distribution = []  # will hold lists of number which fulfill criteria
distribution_carb = []  # will hold lists of number which fulfill criteria

# Need a loop for each distribution in case there are different dimensions
# For Sodium:
for i in range(0, len(distribution_pts)):  # look at each range defined in distribution_pts
    values = []  # will hold all values which meet the criteria on each loop
    for j in range(0, len(dist_one_list)):  # looks through all the sodium ions
        if distribution_pts[i] <= dist_one_list[j] < distribution_pts[i + 1]:  # checks to see if the distance from
            # the pore center is between the two distances defined by distribution_pts
            values.append(dist_one_list[j])  # if the criteria is met, it is added to the list of values
    distribution.append(values)  # the list of values is added to the distribution list as a separate list

# For Carbonyl:
for i in range(0, len(distribution_pts_carb)):  # look at each range defined in distribution_pts
    values = []  # will hold all values which meet the criteria on each loop
    for j in range(0, len(dist_one_list_carb)):  # looks through all the sodium ions
        if distribution_pts_carb[i] <= dist_one_list_carb[j] < distribution_pts_carb[i + 1]:  # checks to see if the distance from
            # the pore center is between the two distances defined by distribution_pts
            values.append(dist_one_list_carb[j])  # if the criteria is met, it is added to the list of values
    distribution_carb.append(values)  # the list of values is added to the distribution list as a separate list

# For Sodium:
group_count = np.zeros((len(distribution_pts), 1))  # will hold the total number of values populating each list in distribution
for i in range(0, len(distribution_pts)):
    total = len(distribution[i])
    group_count[i] = total

# For Carbonyl:
group_count_carb = np.zeros((len(distribution_pts_carb), 1))  # will hold the total number of values populating each list in distribution
for i in range(0, len(distribution_pts_carb)):
    total = len(distribution_carb[i])
    group_count_carb[i] = total

plt.figure(4)
plt.bar(distribution_pts, group_count, bar_width)
#plt.bar(distribution_pts_carb, group_count_carb, bar_width)
plt.title('Distribution of Distances from Center of Pore')
plt.xlabel('Distance from pore center (nm)')
plt.ylabel('Number of ions at this distance')

plt.figure(5)
plt.bar(distribution_pts_carb, group_count_carb, bar_width)
plt.title('Distribution of Distances from Center of Pore')
plt.xlabel('Distance from pore center (nm)')
plt.ylabel('Number of Carbons at this distance')
axes = plt.gca()
axes.set_xlim([0,1.5])

# Density of ions per area

# For Sodium:
areas = np.zeros((len(distribution_pts), 1))
for i in range(0, len(distribution_pts) - 1):
    areas[i] = math.pi*distribution_pts[i + 1]**2 - math.pi*distribution_pts[i]**2

density = np.zeros((len(areas), 1))
for i in range(0, len(areas)):
    if areas[i] == 0:
        density[i] = 0
    else:
        density[i] = group_count[i]/areas[i]
bar_width_2 = bar_width / math.pi
plt.figure(6)
plt.bar(areas, density, bar_width_2)
plt.title('Density of Ions per Area')
plt.xlabel('Area (nm^2) of annulus from r = 0 outwards')
plt.ylabel('Density (Ions/nm^2)')

# For Carbonyl:
areas_carb = np.zeros((len(distribution_pts_carb), 1))
for i in range(0, len(distribution_pts_carb) - 1):
    areas_carb[i] = math.pi*distribution_pts_carb[i + 1]**2 - math.pi*distribution_pts_carb[i]**2

density_carb = np.zeros((len(areas_carb), 1))
for i in range(0, len(areas_carb)):
    if areas_carb[i] == 0:
        density_carb[i] = 0
    else:
        density_carb[i] = group_count_carb[i]/areas_carb[i]
bar_width_3 = bar_width / math.pi
plt.figure(7)
plt.bar(areas_carb, density_carb, bar_width_3)
plt.title('Density of C6 per Area')
plt.xlabel('Area (nm^2) of annulus from r = 0 outwards')
plt.ylabel('Density (Ions/nm^2)')
axes = plt.gca()
axes.set_xlim([0,0.35])



# print 'The average distance of ions from the pore center are as follows:'
# print 'Pore 1: {:2.2f} nm' .format(float(avg[0]))
# print 'Pore 2: {:2.2f} nm' .format(float(avg[1]))
# print 'Pore 3: {:2.2f} nm' .format(float(avg[2]))
# print 'Pore 4: {:2.2f} nm' .format(float(avg[3]))
# print ''
# print 'The standard deviation in distance of ions from the center of the pore are as follows:'
# print 'Pore 1: {:2.2f} nm' .format(float(std[0]))
# print 'Pore 2: {:2.2f} nm' .format(float(std[1]))
# print 'Pore 3: {:2.2f} nm' .format(float(std[2]))
# print 'Pore 4: {:2.2f} nm' .format(float(std[3]))
# print ''
# print 'The average distance of ions from the pore center is %s nm' %Average_Distance
# print 'The average deviation from the pore center is %s nm' %Average_Deviation
# print ''
# print 'Pore Orientation for reference: '
# print '      Pore 1 ------------- Pore 2'
# print '           /            /  '
# print '          /            /  '
# print '         /            /  '
# print '        /            /  '
# print ' Pore 3 ------------- Pore 4'
# print ''
#
#plt.show()

Pore_list = [pore13, pore14, pore23, pore34, pore12, pore42]
Sum_list = []
for i in range(0, len(Pore_list)):
    Sum_list.append([])
    Sum_list[i] = 0


for i in range(0, len(Sum_list)):
    for k in range(traj_start, len(pore12)):
        Sum_list[i] += Pore_list[i][k]

frames = float(1 / float(traj_points - traj_start))
Averages = np.dot(Sum_list, frames).tolist()
Averages.sort()
print Averages
print np.mean(Averages[0:4])

x = x_axis
y = y_axis

pts, traj_pts = np.shape(x)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(-2, 6), ylim=(-2, 6))
line, = ax.plot([], [], lw=2)  # Can include data here

# Create an initialization function in order to hide any plot elements that are not wanted to be shown in every frame
def init():
    line.set_data([], [])
    return line,

# Animate function. This is called sequentially. The purpose is to update the plot elements for each frame
def animate(i):
    x_traj = [x[0][i], x[1][i], x[3][i], x[2][i], x[0][i]]  # [Pore1, Pore2, Pore4, Pore3]
    y_traj = [y[0][i], y[1][i], y[3][i], y[2][i], y[0][i]]
    dist = []
    dist.append(math.sqrt((x_traj[0] - x_traj[1])**2 + (y_traj[0] - y_traj[1])**2))  # Pore1 to Pore2
    dist.append(math.sqrt((x_traj[1] - x_traj[2])**2 + (y_traj[1] - y_traj[2])**2))  # Pore2 to Pore4
    dist.append(math.sqrt((x_traj[3] - x_traj[2])**2 + (y_traj[3] - y_traj[2])**2))  # Pore4 to Pore3
    dist.append(math.sqrt((x_traj[0] - x_traj[3])**2 + (y_traj[0] - y_traj[3])**2))  # Pore1 to Pore3
    line.set_data(x_traj, y_traj)
    annotate1 = ax.annotate(('%1.3f A' %dist[0]), xy=((x_traj[0] + x_traj[1])/2, (y_traj[0] + y_traj[1])/2 + .25))
    annotate2 = ax.annotate(('%1.3f A' %dist[1]), xy=((x_traj[2] + x_traj[1])/2 - 1.25, (y_traj[2] + y_traj[1])/2))
    annotate3 = ax.annotate(('%1.3f A' %dist[2]), xy=((x_traj[3] + x_traj[2])/2 - .25, (y_traj[3] + y_traj[2])/2 - .5))
    annotate4 = ax.annotate(('%1.3f A' %dist[3]), xy=((x_traj[0] + x_traj[2])/2 + 2.25, (y_traj[0] + y_traj[2])/2))
    return line, annotate1, annotate2, annotate3, annotate4


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=traj_pts, interval=200, blit=True)

plt.title('Pore Center Trajectory')
plt.ylabel('Y Position of Pore Center')
plt.xlabel('X Position of each Pore Center')
plt.show()