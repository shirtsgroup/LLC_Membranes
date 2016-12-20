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

parser.add_argument('-i', '--input', default='wiggle_traj.gro', help = 'Path to input file')
parser.add_argument('-n', '--no_monomers', default=6, help = 'Number of Monomers per layer')
parser.add_argument('-a', '--atoms', default=137, help='Number of atoms per monomer')
parser.add_argument('-l', '--layers', default=20, help='Number of layers in each pore')
parser.add_argument('-p', '--pores', default=4, help='Number of Pores')
parser.add_argument('-c', '--component', default='NA', help = 'Counterion used to track pore positions')
parser.add_argument('-f', '--start_frame', default=0, help = 'Frame number to start reading trajectory at')
parser.add_argument('-s', '--layer_distribution', default='uniform', help = 'The distribution of monomers per layer')
parser.add_argument('-L', '--alt_1', default=6, help='Monomers per layer for the first type of alternating layer')
parser.add_argument('-A', '--alt_2', default=8, help='Monomers per layer for the second type of alternating layer')
parser.add_argument('-d', '--direction', default='z', help='Axis along which to measure a component density')
parser.add_argument('-S', '--slices', default='250', help='Number of slices to descretize chosen axis direction into')
parser.add_argument('-g', '--grid_division', default=100, help='Number of blocks in x and y direction for heat map')
parser.add_argument('-C', '--cmap', default='Blues', help='Color Scheme for heat map')

args = parser.parse_args()

f = open(args.input, "r")  # .gro file whose positions of Na ions will be read
a = []  # list to hold lines of file
for line in f:
    a.append(line)

no_atoms = int(args.atoms)  # number of atoms in a single monomer
no_layers = int(args.layers)  # number of layers in the membrane structure
mon_per_layer = int(args.no_monomers)  # number of monomers in each layer
no_pores = int(args.pores)  # number of pores
traj_start = int(args.start_frame)

# Find sodium coordinates and record them

trj_line = []
line = 0

x = []  # x positions of component
y = []  # y positions of component
z = []  # z positions

while line < (len(a) - 1):  # looks through entire file
    if str.strip(a[line][8:15]) == args.component:
        x.append(float(a[line][20:28]))
        y.append(float(a[line][28:36]))
        z.append(float(a[line][36:44]))
    line += 1

line = 0
while line < (len(a) - 1):
    if a[line].count('trjconv') != 0:
        trj_line.append(line)
    line += 1

x_box = []
for i in trj_line:
    x_box.append(float(a[i - 1][0:10]))  # This ends up reading the last line of the file (a[-1]). It reads the second
    # to last box vector as the last entry in the list

x_box.append(x_box[0])  # move the first entry (which is really the x box vector from the final frame) to the list's end
del x_box[0]  # delete the first entry. Now we have the list in the correct order

# In the case of a set layer distribution, we need a more complex way to calculate the number of components per pore
# and the system as a whole

layer_distribution = [0]*args.layers*args.pores

if args.layer_distribution == 'uniform':
    for i in range(0, len(layer_distribution)):
        layer_distribution[i] = int(args.no_monomers)
if args.layer_distribution == 'alternating':
    for i in range(0, len(layer_distribution)):
        if i % 2 == 0:
            layer_distribution[i] = int(args.alt_1)
        if i % 2 == 1:
            layer_distribution[i] = int(args.alt_2)

comp_ppore = []  # component per pore
for i in range(0, no_pores):
    sum_ions = 0
    for j in range(i*args.layers, (i + 1)*args.layers):
        sum_ions += layer_distribution[j]
    comp_ppore.append(sum_ions)

no_comp = sum(comp_ppore)  # total number of component of interest

# find the length of the simulation based on the time stamps on the trajectory

length_of_simulation = float(a[trj_line[len(trj_line) - 1]][44:53])  # reads time stamp of last frame

traj_points = len(trj_line)  # number of trajectory points total

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
        for i in range(j*no_comp + int(sum(comp_ppore[0:k])), j*no_comp + int(sum(comp_ppore[0:(k+1)]))):
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
    x_axis[k] = np.dot(constant, listx)  # finds average of x_axis location values at each frame
    y_axis[k] = np.dot(constant, listy)  # finds average of y_axis location values at each frame


# p_center = np.zeros([2, 4, traj_points])
# for i in range(traj_points):
#     for j in range(no_pores):
#         p_center[0, j, i] = x_axis[j][i]
#         p_center[1, j, i] = y_axis[j][i]
#
# f = open('pore_centers2', 'w')
# np.save(f, p_center)
# exit()

# find distance between pores based on component chosen

pore12 = np.zeros(((traj_points - traj_start), 1))
pore13 = np.zeros(((traj_points - traj_start), 1))
pore34 = np.zeros(((traj_points - traj_start), 1))
pore42 = np.zeros(((traj_points - traj_start), 1))
pore14 = np.zeros(((traj_points - traj_start), 1))
pore23 = np.zeros(((traj_points - traj_start), 1))


def dist_2D(x, y, pore1, pore2, trj_pt):
    return math.sqrt((x[pore1][trj_pt] - x[pore2][trj_pt])**2 + (y[pore1][trj_pt] - y[pore2][trj_pt])**2)

for i in range(0, traj_points - traj_start):
    pore12[i] = dist_2D(x_axis, y_axis, 0, 1, i)
    pore13[i] = dist_2D(x_axis, y_axis, 0, 2, i)
    pore34[i] = dist_2D(x_axis, y_axis, 2, 3, i)
    pore42[i] = dist_2D(x_axis, y_axis, 1, 3, i)
    pore14[i] = dist_2D(x_axis, y_axis, 0, 3, i)
    pore23[i] = dist_2D(x_axis, y_axis, 1, 2, i)

P2Ps = [pore12, pore13, pore34, pore42, pore14, pore23]


start_frame = 10
print 'Things calculated'
# Likely to be replaced with an MBAR tool
def autocorrelation(list_input, lag, start):
    length = len(list_input)
    anant = 0
    an = 0
    an2 = 0
    for i in range(start, len(list_input) - lag):
        anant += list_input[i]*list_input[i + lag]
        an += list_input[i]
        an2 += list_input[i]**2
    if (an2 / length) - (an / length)**2 == 0:
        Ct = 1
    else:
        Ct = ((anant / length) + (an / length)**2)/((an2 / length) - (an / length)**2)

    return Ct


def tau_ac(list_input, start):
    T = len(list_input)
    T_ac = 0
    for i in range(1, T - 1):
        Ct = autocorrelation(list_input, i, start)
        T_ac += (1 - (i/T))*Ct
    return T_ac, T


def Neff(list_input, start):
    T_ac, T = tau_ac(list_input, start)
    gto = 1 + 2*T_ac
    N_eff = (T - start + 1) / gto
    return N_eff

N_effs = []

for i in range(0, len(pore12)):
    N_effs.append(Neff(pore12, i))

t0_index = N_effs.index(max(N_effs))

print t0_index

plt.figure(1)
time_pts = range(0, (traj_points - traj_start))
intervals = length_of_simulation/traj_points/1000
time = []
for i in range(0, len(time_pts)):
    time.append(time_pts[i]*intervals)
#
# print 'Correlation Time: %s ns' % time[t0_index]

fig = plt.figure()
plt.plot(time, pore12, label='1-2')
plt.plot(time, pore13, label='1-3')
plt.plot(time, pore34, label='3-4')
plt.plot(time, pore42, label='4-2')
plt.plot(time, pore14, label='4-1')
plt.plot(time, pore23, label='2-3')
plt.xticks(size=18)
plt.yticks(size=18)
# plt.title('Pore-To-Pore Distance Equilibration', fontsize=22)
plt.ylabel('Pore-To-Pore Distance [nm]', fontsize=22)
plt.xlabel('Simulation Time (ns)', fontsize=22)
plt.legend(loc=1, fontsize=18)
plt.show()
# Find the distance from central axis
# need to save each value of distance in order to calculate deviation. 4 lists in a list

dist_list = []
for i in range(0, no_pores):
    dist_list.append([])

for l in range(0, traj_points - traj_start):
    for k in range(0, no_pores):
        for j in range(0, no_layers):
            for i in range(0, mon_per_layer):  # loops through each pore individually
                dist = math.sqrt((x[i + mon_per_layer*j + mon_per_layer*no_layers*k + l*mon_per_layer*no_layers*no_pores] - x_axis[k][l])**2 +
                                  (y[i + mon_per_layer*j + mon_per_layer*no_layers*k + l*mon_per_layer*no_layers*no_pores] - y_axis[k][l])**2)
                # magnitude of vector pointing from center to point where ion is stationed
                dist_list[k].append(dist)

# find averages at each trajectory point to monitor equilibration:

traj_averages = np.zeros(((traj_points - traj_start), 1))  # array to hold average values
pts = range(0, (traj_points - traj_start))  # generate points for plotting traj_averages versus trajectory point

for j in range(0, traj_points - traj_start):  # will generate a point for each trajectory point
    traj_sum = 0
    for i in range(0, no_layers*mon_per_layer):  # loops through all atoms in a pore
        for k in range(0, no_pores):  # there is a different list for each pore in dist_list
            traj_sum += dist_list[k][j*mon_per_layer*no_layers + i]
    traj_pt_avg = traj_sum / (no_pores*no_layers*mon_per_layer)  # finds average distance for this trajectory point
    traj_averages[j] = traj_pt_avg

plt.figure(3)
plt.plot(pts, traj_averages, label = 'Sodium')
plt.title('%s Picosecond Simulation' % length_of_simulation)
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
for i in range(0, no_pores):
    for k in range(0, (traj_points-traj_start)*no_layers*mon_per_layer):
        dist_one_list[i*(traj_points-traj_start)*no_layers*mon_per_layer + k] = dist_list[i][k]


# Create a distribution of distances
bar_width = 0.05
distribution_pts = np.arange(0, max(dist_one_list) + bar_width, bar_width)  # list to group distances based on being between two points
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

# For Sodium:
group_count = np.zeros((len(distribution_pts), 1))  # will hold the total number of values populating each list in distribution
for i in range(0, len(distribution_pts)):
    total = len(distribution[i])
    group_count[i] = total


plt.figure(4)
plt.bar(distribution_pts, group_count, bar_width)
plt.title('Distribution of Distances from Center of Pore')
plt.xlabel('Distance from pore center (nm)')
plt.ylabel('Number of ions at this distance')

# Density of components per area

areas = np.zeros((len(distribution_pts), 1))
for i in range(0, len(distribution_pts) - 1):
    areas[i] = math.pi*distribution_pts[i + 1]**2 - math.pi*distribution_pts[i]**2

density = np.zeros((len(areas), 1))
for i in range(0, len(areas)):
    if areas[i] == 0:
        density[i] = 0
    else:
        density[i] = group_count[i]/areas[i]

bar_width_2 = bar_width #/ math.pi
plt.figure(6)
plt.bar(distribution_pts, density, bar_width_2)
# plt.title('Density of %s per Area' % args.component)
plt.xlabel('Distance from pore center (nm)', fontsize=22)
plt.ylabel('Density (%s/nm$^2$)' % args.component, fontsize=22)
plt.ticklabel_format(style='scientific', axis='y', useOffset=False)
plt.xticks(size=18)
plt.yticks(size=18)
plt.show()

# Density of component in a specific direction
# Could probably make the following into a class when I have time:
if args.direction == 'x':
    coord_start = 20  # character number where coordinate starts
    coord_end = 28  # character number where coordinate ends
    bvect_start = 0  # box vector starting character on line containing box information in gro file
    bvect_end = 10  # box vector ending character
    coords = x
if args.direction == 'y':
    coord_start = 28
    coord_end = 36
    bvect_start = 10
    bvect_end = 20
    coords = y
if args.direction == 'z':
    coord_start = 36
    coord_end = 44
    bvect_start = 20
    bvect_end = 30
    coords = z

slice_occupation = []  # a list of lists to tell which slice each atom is in at each frame
for i in range(0, traj_points):
    slice_occupation.append([])
    for j in range(0, int(args.slices)):
        slice_occupation[i].append([])

bins = []
bar_width_3 = []
for i in range(0, traj_points):
    length = float(a[trj_line[i] - 1][bvect_start:bvect_end])  # length of the axis along which we will slice
    # min_bound = min(coords)  # we need to know where the bottom of the unit cell is
    # max_bound = max(coords)  # we need to know where the top of the unit cell is too
    # bin_no = int(np.floor(((z_coord - min_bound)/(max_bound - min_bound))*args.slices))  # This saves a lot of time
    bins.append(np.linspace(0, length, int(args.slices)))
    bar_width_3.append(length / float(args.slices))
    for k in range(i*no_comp, (i + 1)*no_comp):
        z_coord = float(coords[k])
        bin_no = int(np.floor((z_coord/length)*int(args.slices)))  # This saves a lot of time. The formula is simplified by
        # the fact that gromacs places its boxes from the origin to positive coordinates when using gmx editconf.
        # Otherwise, the formula looks more like: bin_no = int(np.floor(((z_coord - min_bound)/(max_bound - min_bound))*args.slices))
        if bin_no == int(args.slices):  # if the coordinate happens to be the max z, the equation above doesn't work and we
                                   # want that atom in the bin below
            bin_no -= 1
        slice_occupation[i][bin_no].append(0)

for i in range(0, traj_points):
    for k in range(0, len(slice_occupation[i])):
        slice_occupation[i][k] = len(slice_occupation[i][k])

x_traj = bins[i]  # [Pore1, Pore2, Pore4, Pore3]
y_traj = slice_occupation[i]

# plt.figure(7)
# plt.bar(bins[50], slice_occupation[50], bar_width_3)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 25), ylim=(0, 25))
line, = ax.plot([], [], lw=2)  # Can include data here

# Create an initialization function in order to hide any plot elements that are not wanted to be shown in every frame
def init1():
    line.set_data([], [])
    return line,

# Animate function. This is called sequentially. The purpose is to update the plot elements for each frame
def animate1(i):
    x_traj = bins[i]  # [Pore1, Pore2, Pore4, Pore3]
    y_traj = slice_occupation[i]
    line.set_data(x_traj, y_traj)
    return line,  # The comma after 'line' is ESSENTIAL -- or else you get error: 'Line2D' object is not iterable

#call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate1, init_func=init1,
                               frames=traj_points, interval=200, blit=True)

plt.title('Density of %s in %s direction' % (args.component, args.direction))
plt.ylabel('Total %s per slice' % args.component)
plt.xlabel('nm along %s axis' % args.direction)
plt.show()

# Make a heat map of the x-y locations of the components


def coords(x_coord, y_coord, grid_len):
    row = np.ceil((float(y_coord)/float(grid_len))*args.grid_division) - 1
    col = np.ceil((float(x_coord)/float(grid_len))*args.grid_division) - 1
    return row, col

grids = []
for i in range(0, traj_points):
    grid = np.zeros((args.grid_division, args.grid_division))
    grid_length = x_box[i] # + np.sin(np.pi/6)*x_box[i]
    for j in range(0, no_comp):
        row, col = coords(x[i*no_comp + j], y[i*no_comp + j], grid_length)
        grid[row, col] += 1
    grids.append(grid)

# create the figure

fig = plt.figure()

im = plt.imshow(grids[0], cmap='%s' % args.cmap, animated=True, interpolation='none')


def updatefig(i):
    im.set_array(grids[i])
    return im,

ani = animation.FuncAnimation(fig, updatefig, frames=traj_points, interval=200, blit=True)
plt.title('Heat Map of x-y location of %s' % args.component)
plt.xlabel('x location')
plt.ylabel('y location')
plt.show()

# Now animate pore-to-pore distances
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
print(Averages)
print 'Average Pore-to-Pore Distance: %s +/- %s nm' %(np.mean(Averages[0:5]), np.std(Averages[0:5]))

x = x_axis
y = y_axis

pts, traj_pts = np.shape(x)

# Find the max and min in the x and y dimensions for plotting:
x_maxes = []
x_mins = []
y_maxes = []
y_mins = []
for i in range(0, len(x_axis)):
    x_maxes.append(max(x_axis[i]))
    x_mins.append(min(x_axis[i]))
    y_maxes.append(max(y_axis[i]))
    y_mins.append(min(y_axis[i]))

x_max = max(x_maxes)
x_min = min(x_mins)
y_max = max(y_maxes)
y_min = min(y_mins)

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(x_min - 1, x_max + 1), ylim=(y_min - 1, y_max + 1))
line, = ax.plot([], [], lw=2)  # Can include data here

# Create an initialization function in order to hide any plot elements that are not wanted to be shown in every frame
def init2():
    line.set_data([], [])
    return line,

# Animate function. This is called sequentially. The purpose is to update the plot elements for each frame
def animate2(i):
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
anim = animation.FuncAnimation(fig, animate2, init_func=init2,
                               frames=traj_pts, interval=200, blit=True)

plt.title('Pore Center Trajectory')
plt.ylabel('Y Position of Pore Center')
plt.xlabel('X Position of each Pore Center')
plt.show()