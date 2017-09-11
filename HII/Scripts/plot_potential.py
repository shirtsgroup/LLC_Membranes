#! /usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import Atom_props
import mdtraj as md
import os


def initialize():

    parser = argparse.ArgumentParser(description='Replace charges in .itp file with charges in mol2 file')

    parser.add_argument('-m', '--monomer', default='NAcarb11V', type=str, help='Name of monomer whose charges we are visualizing')
    parser.add_argument('--grid', default=5, help='Number of points in each direction')
    parser.add_argument('-r', '--radius', default=0.1, type=float, help = 'radius of sphere plotted from each point charge')
    parser.add_argument('--atoms', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Which atoms to plot potential of')

    args = parser.parse_args()

    return args


def potential(q, r):
    # V(r) = q / (4 * pi * vacuum permittivity * r)
    return q / (4 * np.pi * epsilon * r)


def potential_grid(pts, q, max_r):
    """
    :param pts: number of pts in the x, y and z direction (assumed x=y=z for now)
    :param q: charge of point at origin
    :return: grid of potential around
    """

    step = max_r / float(pts)

    grid = np.zeros([pts, pts, pts])
    for i in range(pts):
        for j in range(pts):
            for k in range(pts):
                if i == 0 and j == 0 and k == 0:
                    grid[i, j, k] = 0
                else:
                    dist = np.linalg.norm([step * i, j * step, k * step])
                    if dist <= max_r:
                        grid[i, j, k] = potential(q, dist)

    return grid


def xyz_convert(grid, max_r, step):

    x, y, z = grid.shape

    new_grid = np.zeros([x*y*z*8, 3])
    intensity = np.zeros([x*y*z*8])

    shift = [[1, 1, 1], [1, 1, -1], [1, -1, -1], [-1, -1, -1], [-1, -1, 1], [-1, 1, 1], [-1, 1, -1], [1, -1, 1]]
    for image in range(8):
        for i in range(x):
            for j in range(y):
                for k in range(z):
                    if grid[i, j, k] != 0:
                        new_grid[image*x*y*z + i*x*y + j*y + k, :] = [step*i*shift[image][0], step*j*shift[image][1], step*k*shift[image][2]]
                        intensity[image*x*y*z + i*x*y + j*y + k] = grid[i, j, k]

    I = []
    for i in intensity:
        if i != 0:
            I.append(i)

    intensity = np.array(I)

    return new_grid[~np.all(new_grid == 0, axis=1)], intensity  # removes all zero lines


def shift_position(grid, new_position):
    """
    Assume that we are shifting from a grid centered at the origin
    :param grid: old grid that needs to be shifted
    :param new_position: the point where the grid needs to be moved to
    :return: a new grid in a new location
    """

    new_grid = np.zeros(grid.shape)

    for i in range(grid.shape[0]):
                new_grid[i, :] = grid[i, :] + np.array(new_position)

    return new_grid

if __name__ == "__main__":

    args = initialize()

    epsilon = 8.8541878176 * 10 ** -12  # F/m

    all_charges = Atom_props.charges(args.monomer)  # get the charges of all atoms in the molecule
    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in
    t = md.load('%s/../top/HII_Monomer_Configurations/%s.gro' % (location, args.monomer))  # load .gro file of molecule to get coordinates
    keep = [a.index for a in t.topology.atoms if a.name in args.atoms]  # get atoms indices of only atoms we are interested in
    atoms_of_interest = t.atom_slice(keep)  # trajectory of only atoms we want included in this calculation

    charges = [all_charges[a.name] for a in atoms_of_interest.topology.atoms]  # charges of all atoms of interest
    pos = atoms_of_interest.xyz  # coordinates of all atoms of interest

    max_r = args.radius

    for i in range(pos.shape[1]):

        q = charges[i]

        grid = potential_grid(args.grid, q, max_r)

        new_grid, intensity = xyz_convert(grid, max_r, max_r / float(args.grid))

        shifted = shift_position(new_grid, pos[0, i, :])

        if i == 0:
            all_grids = shifted
            all_intensities = intensity
            print all_intensities.shape
        else:
            all_grids = np.concatenate((all_grids, shifted), axis=0)
            all_intensities = np.concatenate((all_intensities, intensity), axis=0)

    # Now need to add all the charges together


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = cm.hsv(intensity/max(intensity))
    colmap = cm.ScalarMappable(cmap=cm.hsv)
    colmap.set_array(all_intensities)
    cb = fig.colorbar(colmap)
    ax.scatter(all_grids[:, 0], all_grids[:, 1], all_grids[:, 2], c=colors, marker='o')
    plt.show()