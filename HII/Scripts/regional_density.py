#! /usr/bin/env python

import mdtraj as md
import argparse
import Structure_char
import numpy as np
from llclib import physical
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the number density of selected groups along axis')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='wiggle.trr', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop doing calculations')
    parser.add_argument('-bins', default=100, type=int, help='Number of bins to use')

    args = parser.parse_args()

    return args


def duplicate(pos, box):
    """
    Duplicate a set of positions periodically once in the +/- xy directions
    :param pos: xyz positions of a set of coordinates to be duplicated periodically
    :param box: box vectors in mdtraj format (t.unitcell_vectors : [nframes, 3, 3]) for every frame
    :return: Periodically duplicated system
    """

    n = pos.shape[1]  # number atoms in original unit cell

    p = np.zeros([pos.shape[0], n*9, 3])  # will hold periodically duplicated system
    p[:, :n, :] = pos

    # x-direction
    for t in range(pos.shape[0]):
        p[t, n:2*n, :] = pos[t, :, :] + box[t, 0, :]
        p[t, 2*n:3*n, :] = pos[t, :, :] - box[t, 0, :]

    # y-direction
    n *= 3
    for t in range(pos.shape[0]):
        p[t, n:2*n, :] = p[t, :n, :] + box[t, 1, :]
        p[t, 2*n:3*n, :] = p[t, :n, :] - box[t, 1, :]

    return p

if __name__ == "__main__":
    
    args = initialize()  # parse the args
    
    print('Loading trajectory...', end="")
    t = md.load('%s' % args.traj, top='%s' % args.gro)[args.begin:args.end]
    print('done')

    regions = ['Tails', 'Head Groups', 'Sodium']
    colors = ['red', 'green', 'blue']
    nT = t.n_frames
    npores = 4

    r_max = 0

    for i, reg in enumerate(regions):

        print('Calculating number density of %s region' % reg)

        pos = Structure_char.restrict_atoms(t, reg)  # restrict trajectory to region

        p = duplicate(pos, t.unitcell_vectors)  # duplicate things periodically

        p_centers = physical.avg_pore_loc(npores, pos)

        equil = 0
        density, r, bin_width = Structure_char.compdensity(pos, p_centers, equil, t.unitcell_vectors, pores=npores, buffer=0, nbins=args.bins)

        plt.bar(r, density, bin_width, color=colors[i], alpha=0.5, label=reg)

    plt.legend()
    plt.ylabel('Component Number Density (number/nm$^2$)')
    plt.xlabel('Distance from pore center (nm)')
    plt.ylim([0, 1.3])
    plt.tight_layout()
    plt.savefig("regional_density.png")
    plt.show()
