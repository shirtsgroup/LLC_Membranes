#!/usr/bin/python

import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import numdensity
from matplotlib import animation


def initialize():

    parser = argparse.ArgumentParser(description='See where the water is at during equilibration')

    parser.add_argument('-f', '--file', default='wiggle_solv.trr', help='Trajectory file (.xtc or .trr should work)')
    parser.add_argument('-c', '--coord', default='wiggle_solv.gro', help='A coordinate file needed by MD traj')
    parser.add_argument('-r', '--radius', default=.6, type=float, help='Radius of cylinder defining pore (nm)')
    parser.add_argument('-b', '--buffer', default=0.1, type=float, help='Percent into membrane to start calculations')
    parser.add_argument('-s', '--save', help='Save the arrays or not', action="store_true")
    parser.add_argument('-l', '--load', help='If youve already save the arrays, load them for speedup', action="store_true")
    parser.add_argument('-B', '--bin', default=0.1, type=float, help='bin size for calculating density')
    parser.add_argument('-g', '--gif', help='Save output as .gif', action="store_true")
    parser.add_argument('-S', '--smooth_factor', type=int, default=5, help='Record measurements every x frames')

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


def cylinder_region(pos, zmin, zmax, radius, pcenters):
    """
    Find the number of water molecules entering into each pore over time
    :param pos: xyz coordinates of all water molecules
    :param zmin: The z dimension of the bottom of the pore for each frame
    :param zmax: The z dimension of the top of the pore for each frame
    :param radius: The radius of a hypothetical cylinder which we are using to define the pore region
    :param pcenters: The locations of the pore centers. Needed to create
    :return: The number of water molecules sucked into the membrane vs. time
    """

    nT = pos.shape[0]
    atomsppore = pos.shape[1] / pcenters.shape[1]
    water = np.zeros([nT])

    for i in range(nT):
        count = 0
        for k in range(pcenters.shape[1]):
            for j in range(atomsppore):
                x, y, z = pos[i, k*atomsppore + j, :]
                if (x - pcenters[0, k, i])**2 + (y - pcenters[1, k, i])**2 <= radius and zmin[i] <= z <= zmax[i]:
                    count += 1
        water[i] = count

    return water


if __name__ == "__main__":

    args = initialize()

    if args.load:
        zmax = np.load('zmax')
        zmin = np.load('zmin')
        pos = np.load('pos')
        time = np.load('time')
        box = np.load('box')
    else:

        t = md.load('%s' % args.file, top='%s' % args.coord)
        time = t.time
        nT = t.xyz.shape[0]
        natoms = t.xyz.shape[1]
        zmax = np.zeros([nT])
        zmin = np.zeros([nT])
        benz_carbs = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
        resnames = np.zeros([natoms])

        # find the thickness of the membrane at each frame
        for i in range(nT):
            z = []
            for a in t.topology.atoms:
                if a.name in benz_carbs:
                    z.append(t.xyz[i, a.index, 2])
            buff = (max(z) - min(z))*args.buffer
            zmax[i] = max(z) - buff
            zmin[i] = min(z) + buff

        box = t.unitcell_lengths  # get the unit cell lengths
        atoms = ['NA']
        atoms_to_keep = [a.index for a in t.topology.atoms if a.name in atoms or a.name == 'O' and 'HOH' in str(a.residue)]  # SOL is stored as HOH in the traj
        t.restrict_atoms(atoms_to_keep)
        pos = t.xyz

        if args.save:
            f = open('zmax', 'w')
            np.save(f, zmax)
            f.close()
            f = open('zmin', 'w')
            np.save(f, zmin)
            f.close()
            f = open('pos', 'w')
            np.save(f, pos)
            f.close()
            f = open('time', 'w')
            np.save(f, time)
            f.close()
            f = open('box', 'w')
            np.save(f, box)

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
    water_counts = cylinder_region(water, zmin, zmax, args.radius, pcenters)
    x, density = numdensity.density(water, 2, args.bin, box, sum='no', smooth_factor=args.smooth_factor)
    plt.figure(1)
    plt.ylabel('Count of waters in pore regions')
    plt.xlabel('Time (ps)')
    plt.plot(time, water_counts)

    # Watch density of water in the z direction evolve over time
    # Use the last frame to determine reasonable bounds on the x and y axes of the following plot

    end = 0
    begin = 0
    while x[-1, begin] < zmin[-1]:
        begin += 1
    while x[-1, end] < zmax[-1]:
        end += 1
    end += 1

    bar_width = x[-1, 1] - x[-1, 0]

    fig = plt.figure(2)
    ax = plt.axes(xlim=(zmin[-1], zmax[-1]), ylim=(0, max(density[-1, begin:end])))

    rects = plt.bar(x[0, :], density[0, :], bar_width, color='c')

    # annotate = ax.annotate('Time: %s ns' % (time[0]/1000), xy=((zmin[-1] + zmax[-1]) / 2, density[-1, end]*2/3))
    # annotate.set_animated(True)

    def init():
        return rects,

    def animate(i):
        """
        http://stackoverflow.com/questions/34372021/python-matplotlib-animate-bar-and-plot-in-one-picture
        """
        for rect, yi in zip(rects, density[i, :]):
            rect.set_height(yi)
            #annotate = ax.annotate('Time: %s ns' % (time[i]/1000.0), xy=((zmin[-1] + zmax[-1]) / 2, density[-1, end]*2/3))
        return rects  #, annotate

    anim = animation.FuncAnimation(fig, animate, frames=int(nT/args.smooth_factor), interval=100, init_func=init)
    plt.ylabel('Count of waters')
    plt.xlabel('Distance into membrane (nm)')
    plt.title('Density of water along z axis')
    if args.gif:
        anim.save('Water_Density.gif', writer='imagemagick')
    plt.show()

