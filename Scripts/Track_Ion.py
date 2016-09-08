#!/usr/bin/python

# Monitor the sodium ions originally in pores to see if they escape into solution when solvated
# The script is used to help decide on the concentration of sodium ions needed in bulk solution while solvating system
# Makes the assumption that the system has already been equilibrated and the membrane thickness is constant

import argparse
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import math
import scipy.optimize as sci

parser = argparse.ArgumentParser(description = 'Monitor the sodium ions originally in pores to see if they escape into \
                                                solution when solvated')
# User inputs with defaults
parser.add_argument('-t', '--trr', default='wiggle.trr', help='Trajectory file')
parser.add_argument('-s', '--tpr', default='wiggle.tpr', help='Atomic level input file needed for gmx trjconv')
parser.add_argument('-g', '--gro', default='wiggle.gro', help='.gro coordinate file for final frame')
args = parser.parse_args()

# Figure out the membrane thickness
output = subprocess.check_output(["Thickness.py", "-i", "%s" % args.gro])
thickness = float(output[20:26])  # membrane thickness

# Run trjconv and get the trajectory information for sodium ions
p1 = subprocess.Popen(["echo", "3"], stdout=subprocess.PIPE)
p2 = subprocess.Popen(["gmx", "trjconv", "-f", "%s" % args.trr, "-s", "%s" % args.tpr, "-o", "wiggle_traj.gro"], stdin=p1.stdout)

p2.wait()

f = open('wiggle_traj.gro', 'r')

a = []
for line in f:
    a.append(line)

line = 0
t = [float(a[line][a[line].index('t=') + 2:len(a[line]) - 1])]  # plus 2 to avoid the '=' sign

while a[line].count('NA') == 0:
    line += 1

z_values = []
z_box_dim = []
while line < len(a) - 1:
    z_values.append([])
    while a[line].count('NA') != 0:
        z_values[len(z_values) - 1].append(float(a[line][36:44]))  # record z position of each sodium ion
        line += 1
    z_box_dim.append(float(a[line][20:30]))  # Get the dimension of the box in the z dimension
    line += 1
    if line >= len(a):  # make sure we aren't at the end of the file
        break  # if we are at the end of the file then don't do the following code block and break the loop
    t.append(float(a[line][a[line].index('t=') + 2:len(a[line]) - 1]))
    line += 1
    while a[line].count('NA') == 0:  # skips through extraneous lines
        line += 1

out_of_bounds = []
for i in range(0, len(z_values)):
    zbound_top = (z_box_dim[i] / 2) + (thickness / 2)  # upper bound of membrane
    zbound_bottom = (z_box_dim[i] / 2) - (thickness / 2)  # lower bound of membrane
    count = 0
    for j in range(0, len(z_values[i])):
        if z_values[i][j] > zbound_top or z_values[i][j] < zbound_bottom:
            count += 1
    out_of_bounds.append(count)

frames = np.arange(1, len(z_values) + 1).tolist()
time = [x/1000 for x in t]

plt.plot(time, out_of_bounds)
plt.xlabel('Time (ns)')
plt.ylabel('Number of Ions out of bounds')
plt.title('Ions diffusing into bulk vs time')
plt.show()