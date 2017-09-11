#!/usr/bin/env python

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import argparse
import bcc_class
from llclib import transform
from llclib import file_rw

def initialize():

    parser = argparse.ArgumentParser(description='Crosslink LLC structure')  # allow input from user

    parser.add_argument('-i', '--input', default='wiggle_init.gro', help = 'Name of input file')
    parser.add_argument('-b', '--build_mon', default='Dibrpyr14.gro', type=str, help='Name of monomer to build with')
    parser.add_argument('-p', '--phase', default='Ia3d', type=str, help='Which infinite minimal surface to build')
    parser.add_argument('-l', '--low', default=0, type=float, help='Lower bound of x, y and z dimension')
    parser.add_argument('--high', default=2*np.pi, type=float, help='Upper bound of x, y and z dimensions')
    parser.add_argument('-n', '--n', default=50, type=int, help='Number of sections to break grid into when '
                                                                'approximating the chosen implicit function')

    args = parser.parse_args()

    return args


def SchwarzD(x):
    """
    :param x: a vector of coordinates (x1, x2, x3)
    :return: An approximation of the Schwarz D "Diamond" infinite periodic minimal surface
    """

    a = np.sin(x[0])*np.sin(x[1])*np.sin(x[2])
    b = np.sin(x[0])*np.cos(x[1])*np.cos(x[2])
    c = np.cos(x[0])*np.sin(x[1])*np.cos(x[2])
    d = np.cos(x[0])*np.cos(x[1])*np.sin(x[2])

    return a + b + c + d


def gyroid(x):

    a = np.sin(x[0])*np.cos(x[1])
    b = np.sin(x[1])*np.cos(x[2])
    c = np.sin(x[2])*np.cos(x[0])

    return a + b + c


def gridgen(surf, low, high, n):

    # make a cubic grid
    x = np.linspace(low, high, n)
    y = np.linspace(low, high, n)
    z = np.linspace(low, high, n)

    if surf == 'Ia3d' or surf == 'gyroid' or surf == 'ia3d':

        gyro = gyroid([x[:, None, None], y[None, :, None], z[None, None, :]])

        gyro_eval = np.zeros([n**3, 3])

        count_gyro = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -0.05 < gyro[i, j, k] < 0.05:
                        gyro_eval[count_gyro, :] = [x[i], y[j], z[k]]
                        count_gyro += 1

        grid = gyro_eval[:count_gyro, :]

    elif surf == 'Pn3m' or surf == 'pn3m':

        schwarz = SchwarzD([x[:, None, None], y[None, :, None], z[None, None, :]])
        schwarz_eval = np.zeros([n**3, 3])
        count_schwarz = 0
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    if -0.05 < schwarz[i, j, k] < 0.05:
                        schwarz_eval[count_schwarz, :] = [x[i], y[j], z[k]]
                        count_schwarz += 1

        grid = schwarz_eval[:count_schwarz, :]
    else:
        print('The phase you selected is not defined (yet)')
        exit()

    return grid


def gradient(v, surf):
    """
    :param v: vector of x, y, z coordinates
    :param phase: which implicit surface is being used to approximate the structure of this phase
    :return: The gradient vector (which is normal to the surface at x)
    """

    x = v[0]
    y = v[1]
    z = v[2]

    if surf == 'Ia3d' or surf == 'gyroid' or surf == 'ia3d':

        a = np.cos(x)*np.cos(y) - np.sin(x)*np.sin(z)
        b = -np.sin(y)*np.sin(x) + np.cos(y)*np.cos(z)
        c = -np.sin(y)*np.sin(z) + np.cos(z)*np.cos(x)

    elif surf == 'Pn3m' or surf == 'pn3m':

        a = np.cos(x)*np.sin(y)*np.sin(z) + np.cos(x)*np.cos(y)*np.cos(z) - np.sin(x)*np.sin(y)*np.cos(z) - np.sin(x)*np.cos(y)*np.sin(z)
        b = np.sin(x)*np.cos(y)*np.sin(z) - np.sin(x)*np.sin(y)*np.cos(z) + np.cos(x)*np.cos(y)*np.cos(z) - np.cos(x)*np.sin(y)*np.sin(z)
        c = np.sin(x)*np.sin(y)*np.cos(z) - np.sin(x)*np.cos(y)*np.sin(z) - np.cos(x)*np.sin(y)*np.sin(z) + np.cos(x)*np.cos(y)*np.cos(z)

    return np.array([a, b, c])


if __name__ == "__main__":

    args = initialize()

    LC = bcc_class.LC(args.build_mon)

    grid = gridgen(args.phase, args.low, args.high, args.n)

    gradv = np.zeros([grid.shape[0], 6])
    for i in range(grid.shape[0]):
        gradv[i, :3] = grid[i, :]
        gradv[i, 3:] = gradient(grid[i, :], args.phase)

    fig = plt.figure(1)
    ax = fig.add_subplot(111, projection='3d')

    X, Y, Z, U, V, W = zip(*gradv)

    ax.quiver(X, Y, Z, U, V, W, length=0.25)
    plt.show()
    exit()
    natoms = LC.natoms

    bcc = np.zeros([natoms*grid.shape[0], 3])

    for i in range(grid.shape[0]):

        n = gradient(grid[i, :], args.phase)  # normal vector to surface at point grid[i, :]

        R = transform.Rvect2vect(LC.linevector, n)  # rotation matrix to rotate monomer in same direction as n
        xyz = transform.rotate_coords(LC.xyz, R)  # rotate all points in bcc monomer with rotation matrix

        # find avg location of reference atoms after rotation (since it will change)
        ref = np.zeros([3])
        for j in range(len(LC.ref_index)):
            ref += xyz[LC.ref_index[j], :]
        ref /= len(LC.ref_index)

        xyz = transform.translate(xyz, ref, grid[i, :])

        bcc[i*natoms:(i + 1)*natoms, :] = xyz

    file_rw.write_gro_pos(bcc, 'initial.gro', res=LC.resid*grid.shape[0], ids=LC.names*grid.shape[0])