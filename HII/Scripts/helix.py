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

layers = 1000
name = 'NA'
radius = 0.2
ionspturn = 7
dbwl = .37
sigma = 0.01

rise = dbwl / ionspturn
print 'rise: %s' % rise

pt = np.array([radius, radius, 0])
pos = np.zeros([layers*ionspturn, 3])

for i in range(layers):
    for j in range(ionspturn):
        Rx = rotate(j*math.pi/(ionspturn/2))
        pos[i*ionspturn + j, :] = np.dot(Rx, pt) + [sigma * np.random.randn(), sigma * np.random.randn(), 0]
        pos[i*ionspturn + j, 2] += (i*ionspturn + j)*rise + sigma * np.random.randn()

f = open('helix.gro', 'w')
f.write('This is a .gro file\n')
f.write('%s\n' % (layers * ionspturn))
for i in range(layers*ionspturn):
    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(i+1, name, name, i+1, pos[i, 0], pos[i, 1], pos[i, 2]) + "\n")

f.write('{:>10f}{:>10f}{:>11f}\n'.format(2*(radius + sigma), 2*(radius + sigma), layers*dbwl + 2*sigma))
f.close()