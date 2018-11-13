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

layers = 20
name = 'NA'
radius = 0.1
ionspturn = 6
dbwl = .37
sigma = 0.01
pos = np.zeros([layers*ionspturn, 3])
pt = np.array([radius, radius, 0])
atoms = layers * ionspturn

for i in range(layers):
    for j in range(ionspturn):
        if i % 2 == 0:
            Rx = rotate(j*math.pi/(ionspturn/2))
        else:
            Rx = rotate((j*math.pi/(ionspturn/2)) + math.pi/ionspturn)
        # Rx = rotate(j*math.pi/(ionspturn/2))
        pos[i*ionspturn + j, :] = np.dot(Rx, pt) + [sigma * np.random.randn(), sigma * np.random.randn(), 0]
        pos[i*ionspturn + j, 2] += i*dbwl + sigma * np.random.randn()

f = open('offset.gro', 'w')
f.write('This is a .gro file\n')
f.write('%s\n' % atoms)
for i in range(atoms):
    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(i+1, name, name, i+1, pos[i, 0], pos[i, 1], pos[i, 2]) + "\n")

f.write('{:>10f}{:>10f}{:>10f}\n'.format(2*(radius + sigma), 2*(radius + sigma), layers*dbwl + 3*sigma))
f.close()
