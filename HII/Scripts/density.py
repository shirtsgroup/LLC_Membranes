#!/usr/bin/env python

import argparse
import mdtraj as md
from Atom_props import mass
import numpy as np
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Measure density of '
                                                 'atoms along a specified axis')
    # User inputs with defaults
    parser.add_argument('-t', '--traj', default='wiggle.trr', help='Trajectory file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help='.gro coordinate file for final frame')
    parser.add_argument('-a', '--atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'])
    parser.add_argument('-A', '--axis', default='z', help='Direction in which to calculate density')
    parser.add_argument('-b', '--bins', default=100, type=int, help='Number of bins for histogram')
    parser.add_argument('-o', '--out', default='density.png', type=str, help='Name of output plot')

    args = parser.parse_args()

    return args


class Density(object):

    def __init__(self, gro, traj, atoms, start=0, end=-1):

        t = md.load(traj, top=gro)[start:end]  # as it stands, this doesn't load last frame

        keep = [a.index for a in t.topology.atoms if a.name in atoms]
        self.mass = [mass[a.name] for a in t.topology.atoms if a.name in atoms]
        self.positions = t.xyz[:, keep, :]
        self.box = t.unitcell_vectors
        self.bin_centers = None
        self.dens = None

        # can implement something here that ensures atoms are in box. For now, just preprocess trajectory with gromacs

    def axial_density(self, axis, bins=100, npores=4):

        axes = {'x': 0, 'y': 1, 'z': 2}
        upper_limit = np.amin(self.box[:, 2, 2]) / 2
        edges = np.linspace(0, upper_limit, bins + 1)
        bin_width = edges[1] - edges[0]  # nm
        self.bin_centers = [(edges[i] + edges[i + 1])/2 for i in range(edges.size - 1)]

        # measure over whole system
        # self.dens = np.zeros([self.positions.shape[0], bins])
        # for t in range(self.positions.shape[0]):
        #     self.dens[t] = np.histogram(self.positions[t, :, axes[axis]], bins=bins, weights=self.mass,
        #                                 range=(0, upper_limit))[0]

        # measure per pore
        atom_per_pore = self.positions.shape[1] // npores
        self.dens = np.zeros([atom_per_pore, npores, bins])
        for t in range(self.positions.shape[0]):
            for p in range(npores):
                self.dens[t, p, :] = np.histogram(self.positions[t, p*atom_per_pore:(p+1)*atom_per_pore, axes[axis]],
                                                  bins=bins, weights=self.mass[p*atom_per_pore:(p+1)*atom_per_pore],
                                                  range=(0, upper_limit))[0]

        self.dens /= bin_width  # units of self.dens are g/mol/nm

        np.savez_compressed('density', dens=self.dens, bins=self.bin_centers)

        return self.dens

    def plot_axial_density(self, save=False, out='density.png'):

        plt.figure()
        data = np.mean(self.dens[:, 2, :], axis=0)
        # data = np.mean(data, axis=0)
        # data /= np.amax(data)
        plt.plot(self.bin_centers, data, linewidth=2)
        plt.ylabel('Number density (grams/mol/nm)', fontsize=14)
        plt.xlabel('Distance (nm)', fontsize=14)
        plt.tight_layout()
        if save:
            plt.savefig(out)
        plt.show()


if __name__ == "__main__":

    args = initialize()

    d = Density(args.gro, args.traj, args.atoms)
    density = d.axial_density(args.axis, bins=args.bins)
    d.plot_axial_density(save=True, out=args.out)