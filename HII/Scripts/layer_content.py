#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from llclib import physical
from matplotlib import animation


def initialize():

    parser = argparse.ArgumentParser(description='Bin atoms into predefined layers to check for uniformity')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Name of trajectory file (.xtc or .trr)')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop calculations')
    parser.add_argument('--skip', default=1, type=int, help='Usage: --skip n . Sample every nth frame')
    parser.add_argument('-a', '--atom', type=str, default='C', help='Name of atoms to track')
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of layers the system was built with')
    parser.add_argument('--shift', default=0, type=float, help='Shift layer centers')

    args = parser.parse_args()

    return args


def assign_pore(pt, p_centers):
    """
    :param pt: xyz coordinates where atom is located
    :param p_centers: coordinates of pore centers (average locations at this frame)
    :param frame: frame
    :return: pore which the point belongs to
    """

    point = np.zeros([2, 1])
    point[:, 0] = pt[:2]

    dist = np.linalg.norm(p_centers - point, axis=0)

    return np.argmin(dist)


def assign_layer(pt, L):
    """
    :param pt: xyz coordinates of point which we want to assign to a layer
    :param L: length of z box vector
    :return: bin
    """

    if pt > L:
        pt -= L

    pt -= args.shift

    return int(np.floor((pt * args.layers) / L))


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[args.begin:args.end:args.skip]
    nT = t.n_frames

    box = t.unitcell_vectors

    keep = [a.index for a in t.topology.atoms if a.name == args.atom]

    pos = t.atom_slice(keep).xyz  # z positions only

    pcenters = physical.avg_pore_loc(4, pos)

    bins = np.zeros([nT, 4*args.layers])
    edges = np.linspace(1, 4*args.layers, 4*args.layers)

    for t in range(nT):
        L = box[t, 2, 2]
        for i in range(pos.shape[1]):
            p = assign_pore(pos[t, i, :], pcenters[:, :, t])
            l = assign_layer(pos[t, i, 2], L)
            # print p, l, pos[t, i, :]
            bins[t, p*args.layers + l] += 1

    # Animated bar plot

    bar_width = 1

    fig = plt.figure(2)
    ax = plt.axes(xlim=(0, 80), ylim=(0, np.mean(bins[0, :]) + 3))

    rects = plt.bar(edges, bins[0, :], bar_width, color='c')

    # annotate = ax.annotate('Time: %s ns' % (time[0]/1000), xy=((zmin[-1] + zmax[-1]) / 2, density[-1, end]*2/3))
    # annotate.set_animated(True)

    def init():
        return rects,

    def animate(i):
        """
        http://stackoverflow.com/questions/34372021/python-matplotlib-animate-bar-and-plot-in-one-picture
        """
        for rect, yi in zip(rects, bins[i, :]):
            rect.set_height(yi)
            #annotate = ax.annotate('Time: %s ns' % (time[i]/1000.0), xy=((zmin[-1] + zmax[-1]) / 2, density[-1, end]*2/3))
        return rects  #, annotate

    anim = animation.FuncAnimation(fig, animate, frames=nT, interval=100, init_func=init)
    plt.ylabel('Count of waters')
    plt.xlabel('Distance into membrane (nm)')
    plt.title('Density of water along z axis')
    # if args.gif:
    #     anim.save('Water_Density.gif', writer='imagemagick')
    plt.show()
    #
    # plt.bar(edges, bins[0, :])
    # plt.show()