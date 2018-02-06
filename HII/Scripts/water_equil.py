#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import numdensity
from matplotlib import animation
from llclib import physical


def initialize():

    parser = argparse.ArgumentParser(description='See where the water is at during equilibration')

    parser.add_argument('-t', '--traj', default='wiggle_solv.trr', help='Trajectory file (.xtc or .trr should work)')
    parser.add_argument('-g', '--gro', default='wiggle_solv.gro', help='A coordinate file needed by MD traj')
    parser.add_argument('-a', '--atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms to based'
                        'membrane thickness on')
    parser.add_argument('-r', '--radius', default=.6, type=float, help='Radius of cylinder defining pore (nm)')
    parser.add_argument('-b', '--buffer', default=0.1, type=float, help='Percent into membrane to start calculations')
    parser.add_argument('-s', '--save', help='Save the arrays or not', action="store_true")
    parser.add_argument('-l', '--load', help='If youve already save the arrays, load them for speedup', action="store_true")
    parser.add_argument('-B', '--bins', default=50, type=int, help='bin size for calculating density')
    parser.add_argument('--gif', help='Save output as .gif', action="store_true")
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

    t = md.load(args.traj, top=args.gro)
    time = t.time
    nT = t.n_frames
    natoms = t.xyz.shape[1]

    keep = [a.index for a in t.topology.atoms if a.name in args.atoms]
    atoms = t.xyz[:, keep, :]
    percent = 5
    n_samples = int(len(keep)*(percent/100))

    # find the thickness of the membrane at last frame
    z = np.sort(atoms[-1, :, 2])
    zmax = np.mean(z[-n_samples:])  # average of the top 'percent' % of z positions in order to suppress outliers
    zmin = np.min(z[:n_samples])
    buff = (zmax - zmin)*args.buffer
    zmax -= buff
    zmin += buff

    water = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']

    hist = np.zeros([nT, args.bins])
    for n in range(nT):
        hist[n, :], edges = np.histogram(t.xyz[n, water, 2], args.bins, range=[zmin, zmax])

    hist /= (zmax - zmin)

    centers = [edges[i] + ((edges[i + 1] - edges[i])/2) for i in range(args.bins)]

    bin_width = centers[1] - centers[0]

    # nT = pos.shape[0]
    # NA = pos[:, :480, :]
    # water = pos[:, 480:, :]
    #
    # filtered_NA = np.zeros([nT, 480, 3])
    # for i in range(nT):
    #     for j in range(480):
    #         if zmin[i] <= NA[i, j, 2] <= zmax[i]:
    #             filtered_NA[i, j, :] = NA[i, j, :]
    #         else:
    #             filtered_NA[i, j, :] = [1000, 1000, 1000]
    #
    # pcenters = avg_pore_loc(4, filtered_NA)
    # water_counts = cylinder_region(water, zmin, zmax, args.radius, pcenters)
    # x, density = numdensity.density(water, 2, args.bin, box, sum='no', smooth_factor=args.smooth_factor)
    # plt.figure(1)
    # plt.ylabel('Count of waters in pore regions')
    # plt.xlabel('Time (ps)')
    # plt.plot(time, water_counts)

    # Watch density of water in the z direction evolve over time
    # Use the last frame to determine reasonable bounds on the x and y axes of the following plot

    fig = plt.figure(2)
    ax = plt.axes(xlim=(zmin, zmax), ylim=(0, 1.1*max(hist[-1, :])))

    rects = plt.bar(centers, hist[0, :], bin_width, color='c')

    def init():
        return rects,

    def animate(i):
        """
        http://stackoverflow.com/questions/34372021/python-matplotlib-animate-bar-and-plot-in-one-picture
        """
        for rect, yi in zip(rects, hist[i, :]):
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

