#!/usr/bin/env python

import mdtraj as md
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib import ticker
import Structure_char


def initialize():

    parser = argparse.ArgumentParser(description='Plot the density of a chosen component as a function of position '
                                                 'within a 2D cross section of a unitcell')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate file')
    parser.add_argument('-a', '--atoms', nargs='+', default=[], help='Reference atom (s) used to locate pores')
    parser.add_argument('-d', '--dim', default='xy', type=str, help='Cross section to look at')
    parser.add_argument('-b', '--bins', default=100, type=int, help='Number of bins')
    parser.add_argument('-r', '--res', nargs='+', default=[], type=str, help='Name of specific residue to look at')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()  # initialize passed arguments

    t = md.load(args.gro)  # load .gro file

    dims = []  # find out how which dimension we are looking at
    for d in args.dim:
        if d == 'x':
            dims.append(0)
        elif d == 'y':
            dims.append(1)
        elif d == 'z':
            dims.append(2)
        else:
            print('Please enter valid dimensions')

    if args.atoms:

        if len(args.atoms) == 1:
            atoms = args.atoms[0]
        else:
            atoms = args.atoms
        # keep only atoms of interest
        pos = Structure_char.restrict_atoms(t, atoms)  # using this function since it has some useful preset groups

    if args.res:

        if len(args.res) == 1:
            res = args.res[0]
            if res == 'SOL':
                res = 'HOH'  # handle mdtraj conversion
        else:
            res = args.res

        keep = [a.index for a in t.topology.atoms if a.residue.name in res]
        pos = t.xyz[:, keep, :]

    H, xedges, yedges = np.histogram2d(pos[0, :, dims[0]], pos[0, :, dims[1]], bins=args.bins)  # bin everything

    # plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    X, Y = np.meshgrid(xedges, yedges)
    if args.atoms:
        if len(args.atoms) == 1:
            plt.title('%s%s density of %s' % (args.dim[0], args.dim[1], atoms))
        else:
            plt.title('%s%s density of %s' % (args.dim[0], args.dim[1], " ".join(atoms)))

    if args.res:
        if len(args.res) == 1:
            plt.title('%s%s density of %s' % (args.dim[0], args.dim[1], res))
        else:
            plt.title('%s%s density of %s' % (args.dim[0], args.dim[1], " ".join(res)))
    plt.xlabel('%s coordinate (nm)' % args.dim[0])
    plt.ylabel('%s coordinate (nm)' % args.dim[1])

    # make a colorbar
    heatmap = ax.pcolormesh(X, Y, H.T, cmap='binary')
    cbar = plt.colorbar(heatmap)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.ax.set_yticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])

    plt.tight_layout()
    plt.show()