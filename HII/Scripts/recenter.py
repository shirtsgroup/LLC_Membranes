#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from llclib import transform
from llclib import file_rw
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Re-position box by rotating or translating it')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of gro file whose box will be modified')
    parser.add_argument('-r', '--rotate', default=7.5, type=float, help='Angle to rotate by')
    parser.add_argument('-t', '--translate', nargs='+', default=[-1, -0.5, 0], type=str, help='[x, y] distances to translate by')
    parser.add_argument('-o', '--output', default='recentered.gro', type=str, help='Name of new .gro produced')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)

    pos = t.xyz[0]
    box = t.unitcell_vectors

    y_height = box[0, 1, 1]
    x_height = box[0, 0, 0]

    # vertices
    v1 = np.zeros(2)
    v2 = box[0][0][:2]
    v3 = box[0][1][:2]
    v4 = v2 + v3

    vertices = np.zeros([4, 3])
    for i, v in enumerate([v1, v2, v3, v4]):
        vertices[i, :2] = v

    avg = [np.mean(vertices[:, 0]), np.mean(vertices[:, 1]), 0]

    t_origin = transform.translate(pos, avg, np.zeros(3))  # translate to the origin

    rotated = transform.rotate_coords_z(t_origin, args.rotate)  # rotate about the origin

    new = transform.translate(rotated, np.zeros(3), avg)  # translate back to where its center was

    move = np.array([float(i) for i in args.translate])

    m = (v3[1] - v1[1]) / (v3[0] - v1[0])
    angle = args.rotate*np.pi /180
    # for i in range(pos.shape[0]):
    #     new[i, :] += move
    #     if new[i, 1] < 0:
    #         new[i, :2] += [y_height*np.sin(angle), y_height*np.cos(angle)]
    #     if new[i, 1] < m*new[i, 0]:
    #         new[i, :2] += [x_height*np.cos(angle), x_height*np.sin(angle)]
    #     # # elif new[i, 0] > y_height:
    #     # #     new[i, 0] -= y_height

    for i in range(pos.shape[0]):
        if new[i, 1] < m*new[i, 0]:
            new[i, :2] += [x_height*np.cos(angle), x_height*np.sin(angle)]

    t.xyz[0] = new

    file_rw.write_gro(t, 'recentered.gro')

    with open(args.gro, 'r') as f:
        a = []
        for line in f:
            a.append(line)