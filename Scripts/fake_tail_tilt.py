#! /usr/bin/env python

import numpy as np
import reposition
import build
import random

# Create benzene ring

benzene = np.zeros([3, 6])
benzene[:, 0] = [-0.576, 0.481, 0.002]
benzene[:, 1] = [-0.553, 0.619, 0.002]
benzene[:, 2] = [-0.660, 0.708, 0.002]
benzene[:, 3] = [-0.791, 0.661, 0.000]
benzene[:, 4] = [-0.815, 0.524, 0.002]
benzene[:, 5] = [-0.708, 0.434, 0.004]

tail_angle = 40  # degrees
tail_angle *= (np.pi / 180)  # converted to radians
length = 11
tail = np.zeros([3, length])
bond_length = .154  # nm

for i in range(1, tail.shape[1]):
    tail[:, i] = tail[:, i - 1] + [bond_length*np.sin(tail_angle), 0, -bond_length*np.cos(tail_angle)]

T = reposition.translate(benzene[:, 0])  # translate to origin

for i in range(benzene.shape[1]):
    cat = np.append(benzene[:, i], [1])
    benzene[:, i] = np.dot(T, cat)[:3]

r = 0.6

T = reposition.translate(benzene[:, 3], r=r)  # translate radially outward

for i in range(benzene.shape[1]):
    cat = np.append(benzene[:, i], [1])
    benzene[:, i] = np.dot(T, cat)[:3]

mpl = 6
dbwl = .37
layers = 20
offset = 'yes'
sigma = 0.5
npores = 4
p2p = 4.5

xyz = np.zeros([3, benzene.shape[1]*mpl*layers*npores])

count = 0
for p in range(npores):  # loop to create multiple pores
    theta = 30  # angle which will be used to do hexagonal packing
    theta*= (np.pi / 180)
    if p == 0:  # unmodified coordinates
        b = 0
        c = 0
    elif p == 1:  # move a pore directly down
        b = -1
        c = 0
    elif p == 2:  # moves pore up and to the right
        b = -np.sin(theta)
        c = -np.cos(theta)
    elif p == 3:  # moves a pore down and to the right
        b = np.cos((np.pi / 2) - theta)
        c = -np.sin((np.pi / 2) - theta)
    for l in range(layers):
        for i in range(mpl):
            if offset == 'yes':
                theta = (np.pi / 3) * i + (l % 2) * (np.pi / 6)
            else:
                theta = (np.pi / 3) * i
            R = reposition.rmatrix(theta)  # rotation matrix to rotate about z axis
            plane = np.array([benzene[:, 0], benzene[:, 2], benzene[:, 4]])  # points defining plane
            T = reposition.translate(benzene[:, 0])  # translate to origin
            Tr = reposition.translate(benzene[:, 3], r=r)  # translate radially outward
            theta = 0
            while theta == 0:
                f = random.choice([-1, 1])  # randomly choose between positive and negative angle
                theta = f * sigma * np.random.randn()
            Rp = build.rotateplane(plane, angle=theta)  # rotate plane in place
            for j in range(6):
                cat = np.append(benzene[:, j], [1])  # translate to origin
                xyz[:, count] = np.dot(T, cat)[:3]
                coord = np.concatenate((xyz[:, count], [1]))
                x = np.dot(Rp, coord)  # rotate plane in place
                xyz[:, count] = x[:3]
                cat = np.append(xyz[:, count], [1])  # translate radially outward
                xyz[:, count] = np.dot(Tr, cat)[:3]
                xyz[:, count] = np.dot(R, xyz[:, count]) + [b*p2p, c*p2p, l*dbwl]

                count += 1

ids = np.array(['C' for id in range(xyz.shape[1])])
res = np.array(['BENZ' for res in range(xyz.shape[1])])

reposition.write_gro(xyz, ids, res, name="benzene_stacked.gro")

xyz = np.zeros([3, tail.shape[1]*mpl*layers*npores])

ids = np.array(['C' for id in range(xyz.shape[1])])
res = np.array(['TAIL' for res in range(xyz.shape[1])])

offset = "no"

count = 0
for p in range(npores):  # loop to create multiple pores
    theta = 30  # angle which will be used to do hexagonal packing
    theta*= (np.pi / 180)
    if p == 0:  # unmodified coordinates
        b = 0
        c = 0
    elif p == 1:  # move a pore directly down
        b = -1
        c = 0
    elif p == 2:  # moves pore up and to the right
        b = -np.sin(theta)
        c = -np.cos(theta)
    elif p == 3:  # moves a pore down and to the right
        b = np.cos((np.pi / 2) - theta)
        c = -np.sin((np.pi / 2) - theta)
    for l in range(layers):
        for i in range(mpl):
            if offset == 'yes':
                theta = (np.pi / 3) * i + (l % 2) * (np.pi / 6)
            else:
                theta = (np.pi / 3) * i
            R = reposition.rmatrix(theta)  # rotation matrix to rotate about z axis
            Tr = reposition.translate(tail[:, 1], r=r)  # translate radially outward
            for j in range(length):
                cat = np.append(tail[:, j], [1])  # translate radially outward
                xyz[:, count] = np.dot(Tr, cat)[:3]
                xyz[:, count] = np.dot(R, xyz[:, count]) + [b*p2p, c*p2p, l*dbwl]

                count += 1

reposition.write_gro(xyz, ids, res, name="tail.gro")


