#! /usr/bin/env python

import mdtraj as md
import argparse
import Structure_char
import numpy as np
from llclib import physical
import matplotlib.pyplot as plt
import os.path as path


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the number density of selected groups along axis')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='wiggle.trr', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop doing calculations')
    parser.add_argument('-bins', default=100, type=int, help='Number of bins to use')
    parser.add_argument('-m', '--multi', nargs='+', help='Overlay the density of each region with the results from '
                                                         'other trajectories')
    parser.add_argument('-s', '--solvate', action="store_true",
                        help='If the system is solvated, plot the number density of water as well')

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

    regions = ['Tails', 'Head Groups', 'Sodium']

    if args.multi:

        colors = ['red', 'green', 'orange', 'blue']

        n = len(args.multi)
        system = np.load(args.multi[0])

        # It is assumed that all of the data uses the same number of bins with the same bin width
        r = system['r']
        bin_width = system['bw']

        results = np.zeros([n, len(regions), len(r)])

        results[0, :, :] = system['results']

        for i in range(1, n):

            system = np.load(args.multi[i])
            results[i, :, :] = system['results']

        for i in range(len(regions)):
            plt.figure(i)
            for j in range(n):
                plt.bar(r, results[j, i, :], bin_width, color=colors[j], alpha=0.6, label=path.splitext(args.multi[j])[0])

            plt.title(regions[i])
            plt.legend()
            plt.ylabel('Component Number Density (number/nm$^2$)')
            plt.xlabel('Distance from pore center (nm)')
            # plt.ylim([0, 0.6])
            plt.tight_layout()
            plt.savefig("%s_density.png" % regions[i])

        plt.show()

        exit()

    colors = ['red', 'green', 'blue']
    print('Loading trajectory...', end="")
    t = md.load('%s' % args.traj, top='%s' % args.gro)[args.begin:args.end]
    print('done')

    nT = t.n_frames
    npores = 4

    r_max = 0

    results = np.zeros([len(regions), args.bins])

    keep = [a.index for a in t.topology.atoms if a.residue.name != 'HOH']  # everything kept if system not solvated

    p_centers = physical.avg_pore_loc(npores, t.xyz[:, keep, :])

    if args.solvate:

        print('Calculating number density of solvent')
        keep = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']
        pos = t.xyz[:, keep, :]
        p = duplicate(pos, t.unitcell_vectors)
        equil = 0
        density, r, bin_width = Structure_char.compdensity(pos, p_centers, equil, t.unitcell_vectors, pores=npores, buffer=0, nbins=args.bins)
        plt.bar(r, density, bin_width, color='orange', alpha=0.5, label='Water')

    for i, reg in enumerate(regions):

        print('Calculating number density of %s region' % reg)

        pos = Structure_char.restrict_atoms(t, reg)  # restrict trajectory to region

        p = duplicate(pos, t.unitcell_vectors)  # duplicate things periodically

        equil = 0
        density, r, bin_width = Structure_char.compdensity(pos, p_centers, equil, t.unitcell_vectors, pores=npores, buffer=0, nbins=args.bins)

        results[i, :] = density

        plt.bar(r, density, bin_width, color=colors[i], alpha=0.5, label=reg)

    np.savez_compressed("density", results=results, r=r, bw=bin_width)

    plt.legend()
    plt.ylabel('Component Number Density (number/nm$^2$)')
    plt.xlabel('Distance from pore center (nm)')
    plt.ylim([0, 1.3])
    plt.tight_layout()
    plt.savefig("regional_density.png")
    plt.show()
