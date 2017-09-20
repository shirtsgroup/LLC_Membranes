#!/usr/bin/env python

import numpy as np
import argparse
import bcc_class
from llclib import transform
from llclib import file_rw
from scipy import spatial
import tqdm
import random
import subprocess
import os


def initialize():

    parser = argparse.ArgumentParser(description='Create bicontinuous cubic membrane structure')  # allow input from user

    parser.add_argument('-b', '--build_mon', default='Dibrpyr14.gro', type=str, help='Name of monomer to build with')
    parser.add_argument('-p', '--phase', default='gyroid', type=str, help='Which infinite minimal surface to build')
    parser.add_argument('-n', '--n', default=50, type=int, help='Number of sections to break grid into when '
                                                                'approximating the chosen implicit function')
    parser.add_argument('-d', '--dim', default=10, type=float, help='Unit cell dimension (length of x, y and z vector)')
    parser.add_argument('-o', '--output', default='initial.gro', type=str, help='Name of output .gro file')
    parser.add_argument('-dens', '--density', default=1.1, type=float, help='Density of system (g/cm3)')
    parser.add_argument('-c', '--curvature', default=-1, type=int,help='-1 : QI phase, 1, QII phase. Determines whether'
                        'the phase is normal or inverted')
    parser.add_argument('-wt', '--weight_percent', default=80.1, type=float, help='Weight %% of monomer in membrane')
    parser.add_argument('-sol', '--solvent', default='glycerol', type=str, help='Name of solvent mixed with monomer')

    args = parser.parse_args()

    return args


def SchwarzD(x):
    """
    :param x: a vector of coordinates (x1, x2, x3)
    :return: An approximation of the Schwarz D "Diamond" infinite periodic minimal surface
    """

    n = 2*np.pi / period

    a = np.sin(n*x[0])*np.sin(n*x[1])*np.sin(n*x[2])
    b = np.sin(n*x[0])*np.cos(n*x[1])*np.cos(n*x[2])
    c = np.cos(n*x[0])*np.sin(n*x[1])*np.cos(n*x[2])
    d = np.cos(n*x[0])*np.cos(n*x[1])*np.sin(n*x[2])

    return a + b + c + d


def gyroid(x):

    n = 2*np.pi / period

    a = np.sin(n*x[0])*np.cos(n*x[1])
    b = np.sin(n*x[1])*np.cos(n*x[2])
    c = np.sin(n*x[2])*np.cos(n*x[0])

    return a + b + c


def gridgen(surf, low, high, n, tol=0.05):

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
                    if -tol < gyro[i, j, k] < tol:
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
                    if -tol < schwarz[i, j, k] < tol:
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

    n = 2*np.pi / period

    if surf == 'Ia3d' or surf == 'gyroid' or surf == 'ia3d':

        a = n*np.cos(n*x)*np.cos(n*y) - n*np.sin(n*x)*np.sin(n*z)
        b = -n*np.sin(n*y)*np.sin(n*x) + n*np.cos(n*y)*np.cos(n*z)
        c = -n*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*z)*np.cos(n*x)

    elif surf == 'Pn3m' or surf == 'pn3m':

        a = n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z)
        b = n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z) - n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z) - n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z)
        c = n*np.sin(n*x)*np.sin(n*y)*np.cos(n*z) - n*np.sin(n*x)*np.cos(n*y)*np.sin(n*z) - n*np.cos(n*x)*np.sin(n*y)*np.sin(n*z) + n*np.cos(n*x)*np.cos(n*y)*np.cos(n*z)

    return np.array([a, b, c])


if __name__ == "__main__":

    args = initialize()

    mw = {'glycerol': 92.09382}

    LC = bcc_class.LC(args.build_mon)  # get all properties of the liquid crystal
    period = args.dim  # length of one side of the unit cell
    grid = gridgen(args.phase, 0, args.dim, args.n)  # 3d grid of points from 0 to args.dim spaced by args.dim / args.n

    # figure out how many grid points we actually want to keep
    NA = 6.022 * 10 ** 23  # avogadros number
    mass = LC.MW / NA  # mass of a single monomer (g)
    V = args.dim ** 3 * 10 ** -21  # volume of unit cell (cm ^ 3)
    mon_mass = (args.weight_percent / 100) * args.density * V  # mass of all monomers in unit cell
    nmon = int(mon_mass / mass)

    mass = mw[args.solvent] / NA
    solv_mass = args.density * V * (100 - args.weight_percent) / 100
    nsol = int(solv_mass / mass)

    while grid.shape[0] > nmon:
        delete = random.randint(0, grid.shape[0] - 1)
        grid = np.delete(grid, delete, 0)

    # print('Filtering grid points to obtain correct density')
    # nn = {}
    # for i in tqdm.tqdm(range(grid.shape[0])):
    #     others = np.delete(grid, i, 0)  # delete self entry or else the nearest neighbor is itself
    #     n = spatial.KDTree(others).query(grid[i, :])[1]  # find index of nearest neighbor
    #     index = np.where(grid == others[n])  # find the index in others corresponding to nn in arr
    #     ld = np.linalg.norm(grid[i, :] - grid[index[0][0]])  # calc linear distance between position and nearest neighbor
    #     nn[i] = ld
    #     if i == 7:
    #         print(n)
    #         print(ld)
    #         exit()

    # nn_sorted = sorted(nn, key=nn.get)[-nmon:]
    # keep = np.zeros([nmon, 3])
    # for i in range(nmon):
    #     keep[i, :] = grid[nn_sorted[i], :]
    #
    # grid = keep

    natoms = LC.natoms  # number of atoms in monomer

    bcc = np.zeros([natoms*grid.shape[0], 3])
    count = 0

    for i in range(grid.shape[0]):

        n = args.curvature * gradient(grid[i, :], args.phase)  # normal vector to surface at point grid[i, :]

        pts = np.zeros([10, 3])
        for k in range(10):
            pts[k, :] = (k/10)*n

        normal = transform.translate(pts, pts[0, :], grid[i, :])

        R = transform.Rvect2vect(LC.linevector, n)  # rotation matrix to rotate monomer in same direction as n

        # translate to origin
        xyz_origin = transform.translate(LC.xyz, LC.reference, np.array([0, 0, 0]))

        xyz_origin = transform.rotate_coords(xyz_origin, R)  # rotate all points in bcc monomer with rotation matrix

        # find avg location of reference atoms after rotation (since it will change)
        ref = np.zeros([3])
        for j in range(len(LC.ref_index)):
            ref += xyz_origin[LC.ref_index[j], :]
        ref /= len(LC.ref_index)

        xyz_origin = transform.translate(xyz_origin, ref, grid[i, :])  # move monomer to grid point w.r.t. reference point on monomer

        bcc[i*natoms:(i + 1)*natoms, :] = xyz_origin

    file_rw.write_gro_pos(bcc, 'initial.gro', res=LC.resid*grid.shape[0], ids=LC.names*grid.shape[0], box=[args.dim, args.dim, args.dim])

    # solvate system

    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    subprocess.call(["gmx", "insert-molecules", "-f", "initial.gro", "-ci", "%s/../top/structures/%s.gro" % (location, args.solvent), "-nmol", "%s" % nsol, "-o", "initial.gro"])