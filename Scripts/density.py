#!/usr/bin/python

import argparse
import subprocess
from Atom_props import mass
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
parser.add_argument('-b', '--buffer', default='1', help='Distance into membrane from min and max where density '
                                                        'measurements will be taken')
args = parser.parse_args()

# Figure out the membrane thickness
output = subprocess.check_output(["Thickness.py", "-i", "%s" % args.gro])
thickness = float(output[20:26])  # membrane thickness
print thickness
#Run trjconv and get the trajectory information for sodium ions
# p1 = subprocess.Popen(["echo", "0"], stdout=subprocess.PIPE)
# p2 = subprocess.Popen(["gmx", "trjconv", "-f", "%s" % args.trr, "-s", "%s" % args.tpr, "-o", "wiggle_traj.gro"], stdin=p1.stdout)

# p2.wait()

f = open('wiggle_traj_bak.gro', 'r')

a = []
for line in f:
    a.append(line)

line = 0
t = [float(a[line][a[line].index('t=') + 2:len(a[line]) - 1])]  # plus 2 to avoid the '=' sign
while a[line].count('HII') == 0:
    line += 1

count = line
res = []  # generate a list of residue names
while a[count].count('Generated') == 0:
    residue = str.strip(a[count][5:8])
    if residue not in res:
        res.append(residue)
    count += 1

del res[len(res) - 1]  # One of the box dimensions gets recorded so delete that

z_limits = []  # top and bottom bound, in between which density will be calculated from coordinates
box_dim = []  # all components of box dimensions at each frame
in_bounds_atoms = []
while line < len(a) - 1:
    z_limits.append([])
    in_bounds_atoms.append([])
    atoms = 0
    while str.strip(a[line][5:8]) in res:  # skip to the end of this frame for now
        line += 1
        atoms += 1
    box_dim.append(a[line])  # Get the dimensions of all components of the box
    z_dim = float(box_dim[len(box_dim) - 1][20:29])  # now pick out just the z-dimension
    upper_lim = (thickness / 2) + (z_dim / 2) - float(args.buffer)  # calculate an upper bound in z for density measurements
    lower_lim = -(thickness / 2) + (z_dim / 2) + float(args.buffer)  # calculate a lower bound in z for density measurements
    z_limits[len(z_limits) - 1] = [upper_lim, lower_lim]  # add these limits to z_limits as we will need them later
    line -= atoms  # go back to the beginning of the frame
    while str.strip(a[line][5:8]) in res:
        if upper_lim >= float(a[line][36:44]) >= lower_lim: #<= upper_lim and float(a[line][36:44]) >= lower_lim:
            in_bounds_atoms[len(in_bounds_atoms) - 1].append(str.strip(a[line][11:15]))  # add qualified atoms' identity
        line += 1
    line += 1
    if line >= len(a):  # make sure we aren't at the end of the file
        break  # if we are at the end of the file then don't do the following code block and break the loop
    t.append(float(a[line][a[line].index('t=') + 2:len(a[line]) - 1])) # get time of trajectory
    line += 1
    while str.strip(a[line][5:8]) not in res:  # skips through extraneous lines
        line += 1

# find the mass of the enclosed system

avagadro = 6.022 * 10 ** 23
tot_mass_per_frame = []
for i in range(0, len(z_limits)):
    tot_mass = 0
    for j in range(0, len(in_bounds_atoms[i])):
        tot_mass += mass[in_bounds_atoms[i][j]]
    tot_mass_per_frame.append(tot_mass / avagadro)

volume = []
for i in range(0, len(box_dim)):
    x = math.sqrt(float(box_dim[i][0:10])**2 + float(box_dim[i][30:40]) ** 2 + float(box_dim[i][40:50]) ** 2)
    y = math.sqrt(float(box_dim[i][10:20])**2 + float(box_dim[i][50:60]) ** 2 + float(box_dim[i][60:70]) ** 2)
    z = abs(z_limits[i][0] - z_limits[i][1])
    volume.append(x*y*z)

nm_to_cm = 10 ** -21
density = []
for i in range(0, len(volume)):
    density.append(tot_mass_per_frame[i] / (volume[i] * nm_to_cm))

t = [x/1000 for x in t]

start = 2* len(density) / 3
print start
print 'Average density = %s' % np.mean(density[start:len(density) - 1])
print 'Standard Deviation Density = %s' % np.std(density[start:len(density) - 1])

plt.plot(t, density)
plt.ylabel('Membrane Density (g/cm$3$)')
plt.xlabel('Time (ns)')
plt.title('Density vs. Time')
plt.show()