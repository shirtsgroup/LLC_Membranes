#!/usr/bin/python

import numpy as np
import math


def rotate(theta):

    Rx = np.zeros((3, 3))  # makes a 3 x 3 zero matrix
    Rx[0, 0] = math.cos(theta)  # This line and subsequent edits to Rx fills in entries needed for rotation matrix
    Rx[1, 0] = math.sin(theta)
    Rx[0, 1] = -math.sin(theta)
    Rx[1, 1] = math.cos(theta)
    Rx[2, 2] = 1

    return Rx

atoms = 600
name = 'NA'
pos = np.zeros([atoms, 3])
radius = 0.1
rise = 0.37
pt = np.array([radius, radius, 0])
layers = 100

for i in range(layers):
    for j in range(6):
        Rx = rotate(j*math.pi/3)
        pos[i*6 + j, :] = np.dot(Rx, pt)
        pos[i*6 + j, 2] += i*rise

f = open('stacked.gro', 'w')
f.write('This is a .gro file\n')
f.write('%s\n' % atoms)
for i in range(atoms):
    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(i+1, name, name, i+1, pos[i, 0], pos[i, 1], pos[i, 2]) + "\n")

f.write('0.000000 0.000000 0.000000\n')
f.close()
