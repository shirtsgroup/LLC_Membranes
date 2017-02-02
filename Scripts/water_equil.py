#!/usr/bin/python

import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='See where the water is at during equilibration')

    parser.add_argument('-f', '--file', default='wiggle_solv.trr', help='Trajectory file (.xtc or .trr should work)')
    parser.add_argument('-c', '--coord', default='wiggle_solv.gro', help='A coordinate file needed by MD traj')
    parser.add_argument('-r', '--radius', default=.5, type=float, help='Radius of cylinder defining pore')
    parser.add_argument('-b', '--buffer', default=0.1, type=float, help='Percent into membrane to start calculations')

    args = parser.parse_args()

    return args


def avg_pore_loc(npores, pos):
    """
    :param no_pores: the number of pores in the unit cell
    :param pos: the coordinates of the component(s) which you are using to locate the pore centers
                      (numpy array with dimensions: [3 (xyz coordinates), no components, no frames])
    :return: numpy array containing the x, y coordinates of the center of each pore at each frame
    """

    # Find the average location of the pores w.r.t. x and y
    nT = pos.shape[0]
    comp_ppore = pos.shape[1] / npores

    p_center = np.zeros([2, npores, nT])

    for i in range(nT):
        for j in range(npores):
            out_of_bounds = 0
            for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                if pos[i, k, 0] == 1000:  # exclusion
                    out_of_bounds += 1
                else:
                    p_center[:, j, i] += pos[i, k, 0:2]
            p_center[:, j, i] /= (comp_ppore - out_of_bounds)  # take the average

    return p_center


def cylinder_region(pos, zmin, zmax, radius):
    """
    Find the number of water molecules entering into each pore over time
    :param pos: xyz coordinates of all water molecules
    :param zmin: The z dimension of the bottom of the pore for each frame
    :param zmax: The z dimension of the top of the pore for each frame
    :param radius: The radius of a hypothetical cylinder which we are using to define the pore region
    :return: The number of water molecules sucked into the membrane vs. time
    """

    nT = pos.shape[0]
    natoms = pos.shape[1]
    water = np.zeros([nT])

    for i in range(nT):
        count = 0
        for j in range(natoms):
            x, y, z = pos[i, j, :]
            if x**2 + y**2 <= radius and zmin[i] <= z <= zmax[i]:
                count += 1
        water[i] = count

    return water

if __name__ == "__main__":
    args = initialize()

    # t = md.load('%s' % args.file, top='%s' % args.coord)
    # nT = t.xyz.shape[0]
    # natoms = t.xyz.shape[1]
    # zmax = np.zeros([nT])
    # zmin = np.zeros([nT])
    # benz_carbs = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
    # resnames = np.zeros([natoms])
    #
    # # find the thickness of the membrane at each frame
    # for i in range(nT):
    #     z = []
    #     for a in t.topology.atoms:
    #         if a.name in benz_carbs:
    #             z.append(t.xyz[i, a.index, 2])
    #     buff = (max(z) - min(z))*args.buffer
    #     zmax[i] = max(z) - buff
    #     zmin[i] = min(z) + buff
    #
    # atoms = ['NA']
    # atoms_to_keep = [a.index for a in t.topology.atoms if a.name in atoms or 'HOH' in str(a.residue)]  # SOL is stored as HOH in the traj
    # t.restrict_atoms(atoms_to_keep)
    # pos = t.xyz
    #
    # f = open('zmax', 'w')
    # np.save(f, zmax)
    # f.close()
    # f = open('zmin', 'w')
    # np.save(f, zmin)
    # f.close()
    # f = open('pos', 'w')
    # np.save(f, pos)
    # f.close()
    zmax = np.load('zmax')
    zmin = np.load('zmin')
    pos = np.load('pos')
    nT = pos.shape[0]
    NA = pos[:, :480, :]
    water = pos[:, 480:, :]

    filtered_NA = np.zeros([nT, 480, 3])
    for i in range(nT):
        for j in range(480):
            if zmin[i] <= NA[i, j, 2] <= zmax[i]:
                filtered_NA[i, j, :] = NA[i, j, :]
            else:
                filtered_NA[i, j, :] = [1000, 1000, 1000]

    pcenters = avg_pore_loc(4, filtered_NA)
    water_counts = cylinder_region(water, zmin, zmax, args.radius)
    t = np.linspace(0, nT - 1, nT)
    plt.plot(t, water_counts)
    plt.show()

