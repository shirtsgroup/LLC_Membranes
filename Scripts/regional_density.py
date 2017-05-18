#! /usr/bin/env python

import mdtraj as md
import argparse
import Structure_char
import numpy as np
from llclib import physical
import copy
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the number density of selected groups along axis')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='wiggle.trr', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-e', '--end', type=int, help='Frame to stop doing calculations')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()  # parse the args

    if args.end:
        t = md.load('%s' % args.traj, top='%s' % args.gro)[args.begin:args.end]
    else:
        t = md.load('%s' % args.traj, top='%s' % args.gro)[args.begin:]

    # Check for certain special arguments
    regions = ['tails', 'benzene', 'NA']
    colors = ['red', 'green', 'blue']
    nT = t.n_frames

    r_max = 0

    for i in range(len(regions)):

        traj = copy.deepcopy(t)  # preserve original trajectory
        pos = Structure_char.restrict_atoms(t, regions[i])  # restrict trajectory to region

        tot_atoms = np.shape(pos)[1]  # number of atoms in the restricted trajectory
        n_pores = 4  # number of pores
        comp_ppore = tot_atoms/n_pores  # number of components in each pore

        p_centers = physical.avg_pore_loc(n_pores, pos)

        equil = 0
        density, r, bin_width = Structure_char.compdensity(pos, p_centers, equil, t.n_atoms, n_pores, buffer=0)
        t = traj

        nonzero = np.trim_zeros(density, trim='b')

        if r[len(nonzero) - 1] > r_max:
            r_max = r[len(nonzero) - 1]

        # plt.bar(r[:stop], density[:stop], bin_width, color=colors[i], alpha=0.5)
        plt.bar(r, density, bin_width, color=colors[i], alpha=0.5, label=regions[i])

        print '%s region calculated' % regions[i]

    plt.legend()
    plt.xlabel('Density')
    plt.ylabel('Distance from pore center (nm)')
    plt.xlim([0, r_max])
    plt.show()