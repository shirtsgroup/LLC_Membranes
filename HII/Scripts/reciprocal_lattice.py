#!/usr/bin/python

import Get_Positions
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import SAXS
import Periodic_Images
import Radial_int_pixels
import mdtraj as md
import random


def hkl(max):

    n = 1  # number of combinations at this max h value
    ng = [n]
    tot = 1
    for i in range(max):
        n += 2 + i
        ng.append(n)
        tot += n

    hkl = np.zeros([tot, 3])

    count = 0
    for i in range(max + 1):
        for j in range(ng[i]):
            hkl[count, 0] = i
            count += 1

    count = 0
    hmax = 1
    for j in ng:
        for i in range(hmax + 1):
            for k in range(i):
                hkl[count, 1] = i - 1
                count += 1
        hmax += 1

    count = 0
    hmax = 0
    for j in ng:
        count += 1
        for i in range(hmax):
            count += 1
            for k in range(i + 1):
                hkl[count, 2] = k + 1
                count += 1
        hmax += 1

    return hkl


def all_hkl(max):

    all = np.zeros([max ** 3, 3])

    count = 0
    for i in range(max):
        for j in range(max):
            for k in range(max):
                all[count, :] = [i, j, k]
                count += 1

    return all


def monoclinic_d(hkl, abc, angle):

    a, b, c = abc
    h, k, l = hkl
    angle *= (np.pi / 180)
    A = h**2 / (a **2 * np.sin(angle)**2)
    B = k**2 / b ** 2
    C = l ** 2 / (c**2 * np.sin(angle)**2)
    D = (2 * h * l * np.cos(angle)) / (a * c * np.sin(angle)**2)

    return np.sqrt(1 / (A + B + C + D))


def cubic_d(hkl, abc, angle):

    h, k, l = hkl

    return abc[0] / np.sqrt(h**2 + k**2 + l**2)


def dhkl(indices, abc, a):

    points = indices.shape[0]

    d = np.zeros([points])

    for i in range(points):
        # d[i] = monoclinic_d(indices[i, :], abc, a)
        d[i] = cubic_d(indices[i, :], abc, a)

    return d


def twotheta(d, l):

    theta = np.zeros(d.shape)

    for i in range(d.shape[0]):
        t = l/(2 * d[i])
        if -1 < t < 1:
            theta[i] = 2*np.arcsin(t) * (180 / np.pi)

    return theta


abc = [9, 9, 4]  # unit cell dimensions [angstroms]
#
# planes = 1
# indices = hkl(planes)
# d = dhkl(indices, abc, 60)
# theta = twotheta(d, 1.54)
# while theta[-1] != 0:
#     planes += 1
#     indices = hkl(planes)
#     d = dhkl(indices, abc, 60)
#     theta = twotheta(d, 1.54)
hlim = (abc[0] / (1.54 / 2))*np.sin(np.pi / 3)

indices = all_hkl(int(hlim))
d = dhkl(indices, abc, 60)
theta = twotheta(d, 1.54)

count = 0
for i in range(theta.shape[0]):
    if theta[i] != 0:
        count += 1
print count
sphere = np.zeros([count, 3])
count = 0
for i in range(theta.shape[0]):
    if theta[i] != 0:
        sphere[count] = indices[i, :]
        count += 1

new_sphere = np.zeros([count*8, 3])
mult = [[1, 1, 1], [1, 1, -1], [1, -1, -1], [-1, -1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [-1, 1, 1]]
for i in range(8):
    for j in range(count):
        new_sphere[count*i + j] = sphere[j, :] * mult[i]

# indices = all_hkl(4)

# theta = twotheta(d, 1.54)
print theta.shape
print 'Maximum d-spacing you can see with this configuration: %s' % max(d[1:])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# ax.scatter(indices[:, 0], indices[:, 1], zs=indices[:, 2])
ax.scatter(new_sphere[:, 0], new_sphere[:, 1], zs=new_sphere[:, 2])

plt.show()