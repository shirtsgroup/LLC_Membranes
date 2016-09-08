#!/usr/bin/python

import argparse
import subprocess
from Atom_props import charge
from Plot_xvg import extract_data
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

start = time.time()

# Calculate the charge distribution across the membrane and in the surrounding water

parser = argparse.ArgumentParser(description = 'Monitor charge across membrane')
# User inputs with defaults
parser.add_argument('-t', '--trr', default='wiggle.trr', help='Trajectory file')
parser.add_argument('-s', '--tpr', default='wiggle.tpr', help='Atomic level input file needed for gmx trjconv')
parser.add_argument('-g', '--gro', default='wiggle.gro', help='.gro coordinate file for final frame')
parser.add_argument('-e', '--edr', default='wiggle.edr', help='Energy file')
parser.add_argument('-d', '--slices', default=100000, help='width of slices to be taken in the z direction when descretizing membrane')

args = parser.parse_args()

p1 = subprocess.Popen(["echo", "0"], stdout=subprocess.PIPE)
p2 = subprocess.Popen(["gmx", "trjconv", "-f", "%s" % args.trr, "-s", "%s" % args.tpr, "-o", "wiggle_traj.gro"], stdin=p1.stdout)
p2.wait()

p3 = subprocess.Popen(["echo", "20"], stdout=subprocess.PIPE)  # 20 gives the z-box dimension
p4 = subprocess.Popen(["gmx", "energy", "-f", "%s" % args.edr], stdin=p3.stdout)
p4.wait()

f = open('wiggle_traj.gro', 'r')

a = []
for line in f:
    a.append(line)

line = 0
t = [float(a[line][a[line].index('t=') + 2:len(a[line]) - 1])]  # plus 2 to avoid the '=' sign

line += 1
no_atoms = int(a[line])
line += 1

count = line
res = []  # generate a list of residue names
while a[count].count('Generated') == 0:
    residue = str.strip(a[count][5:8])
    if residue not in res:
        res.append(residue)
    count += 1

del res[len(res) - 1]  # One of the box dimensions gets recorded so delete that

t = extract_data('energy.xvg')[0]  # Uses Plot_xvg.py to get time points as a list
z_box = extract_data('energy.xvg')[1]  # uses Plot_xvg.py to get z_box dimensions as a list

no_frames = len(t)  # number of trajectory frames

slice_occupation = []  # a list of lists to tell which slice each atom is in at each frame
for i in range(0, no_frames):
    slice_occupation.append([])
    for j in range(0, args.slices):
        slice_occupation[i].append([])

z_min = []
z_max = []
for i in range(0, no_frames):
    z_vals = []  # make a list of all the z_values in this frame
    count = line  # make a new variable to count lines we go through
    while str.strip(a[count][5:8]) in res:  # continue looping until we are no longer reading atom information
        z_vals.append(float(a[count][36:44]))
        count += 1
    min_bound = min(z_vals)  # we need to know where the bottom of the unit cell is
    max_bound = max(z_vals)  # we need to know where the top of the unit cell is too
    z_min.append(min_bound)
    z_max.append(max_bound)
    while str.strip(a[line][5:8]) in res:
        z_coord = float(a[line][36:44])
        bin_no = int(np.floor(((z_coord - min_bound)/(max_bound - min_bound))*args.slices))  # This saves a lot of time
        if bin_no == args.slices:  # if the coordinate happens to be the max z, the equation above doesn't work and we
                                   # want that atom in the bin below
            bin_no -= 1
        # if str.strip(a[line][8:15]) == 'NA':
        #     slice_occupation[i][bin_no].append(charge['NA'])
        slice_occupation[i][bin_no].append(charge[str.strip(a[line][8:15])])
        line += 1
    min_bound = max_bound  # now the lower bound is the old upper bound.. and we proceed upwards through the box
    while str.strip(a[line][5:8]) not in res:
        line += 1
        if line >= len(a):  # make sure we aren't at the end of the file
            break  # if we are at the end of the file then don't do the following code block and break the loop
    for k in range(0, len(slice_occupation[i])):
        slice_occupation[i][k] = sum(slice_occupation[i][k])

end = time.time()

print end - start

z = []
for i in range(0, no_frames):
    z.append([])
    for j in range(0, args.slices):
        z[i].append([])
        z[i][j] = z_min[i] + j*(z_max[i]/args.slices)

plt.plot(z[30], slice_occupation[30])
plt.ylabel('Charge in slice of membrane')
plt.xlabel('Z Box Dimension (nm)')
plt.title('Charge Distribution')

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = plt.axes(xlim=(0, 25), ylim=(-25, 25))
line, = ax.plot([], [], lw=2)  # Can include data here

# Create an initialization function in order to hide any plot elements that are not wanted to be shown in every frame
def init():
    line.set_data([], [])
    return line,

# Animate function. This is called sequentially. The purpose is to update the plot elements for each frame
def animate(i):
    x = np.array(z[i])
    y = np.array(slice_occupation[i])
    # x = np.linspace(0, 2, 1000)
    # y = np.sin(2 * np.pi * (x - 0.01 * i))
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=no_frames, interval=200, blit=True)

plt.title('Pore Center Trajectory')
plt.ylabel('Y Position of Pore Center')
plt.xlabel('X Position of each Pore Center')
plt.show()