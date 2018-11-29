#!/usr/bin/env python

import argparse
import mdtraj as md
from LLC_Membranes.llclib import physical, topology
import numpy as np
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate radial distribution of a monomer component')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='traj_whole.xtc', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-r', '--residue', default='ETH', help='Name of residue whose radial distribution function we'
                                                               'want to calculate')
    parser.add_argument('-b', '--build_monomer_residue', default='HII', help='Name of monomer residue used to build'
                                                                             'HII phase')
    parser.add_argument('-bins', '--bins', default=50, type=int, help='Number of bins for histogram of distances')
    parser.add_argument('--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('--end', default=-1, type=int, help='Frame to stop doing calculations')
    parser.add_argument('-skip', '--skip', default=1, type=int, help='Analyze every skip frames')
    parser.add_argument('-l', '--load', help='Name of compressed .npz to load')
    parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials')

    return parser


class System(object):

    def __init__(self, gro, traj, residue, monomer, begin=0, end=-1, skip=1, npores=4):

        self.t = md.load(traj, top=gro)[begin:end:skip]
        self.npores = npores

        self.residue = topology.Residue(residue)
        self.monomer = topology.LC('%s.gro' % monomer)
        res = [a.index for a in self.t.topology.atoms if a.residue.name == residue]
        self.com = physical.center_of_mass(self.t.xyz[:, res, :], [v for v in self.residue.mass.values()])

        self.r = None
        self.density = None
        self.bootstraps = None
        self.errorbars = None

    def radial_distribution_function(self, bins=50):

        self.r = np.zeros([bins])
        self.density = np.zeros([self.t.n_frames, bins])

        pore_defining_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms]
        pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, pore_defining_atoms, :])
        self.r, self.density = physical.compdensity(self.com, pore_centers, self.t.unitcell_vectors, pores=self.npores,
                                                    nbins=bins, rmax=3.5)

    def bootstrap(self, nboot):

        nT = self.density.shape[0]
        self.bootstraps = np.zeros([nboot, self.density.shape[1]])

        for i in range(nboot):
            frames = np.random.choice(nT, size=nT, replace=True)  # choose random frames with replacement
            self.bootstraps[i, :] = self.density[frames, :].mean(axis=0)  # average rdf from all frames

        confidence = 95  # percent confidence interval
        lower_confidence = (100 - confidence) / 2
        upper_confidence = 100 - lower_confidence

        self.errorbars = np.zeros([2, self.density.shape[1]])
        self.errorbars[0, :] = np.abs(np.percentile(self.bootstraps, lower_confidence, axis=0) -
                                      self.density.mean(axis=0))  # 2.5 percent of data below this value
        self.errorbars[1, :] = np.percentile(self.bootstraps, upper_confidence, axis=0) - self.density.mean(axis=0)

    def plot(self):

        if self.bootstraps is not None:
            plt.plot(self.r, self.bootstraps.mean(axis=0), linewidth=2)
            plt.fill_between(self.r, self.errorbars[1, :] + self.density.mean(axis=0), self.density.mean(axis=0) -
                             self.errorbars[0, :], alpha=0.7)
        else:
            plt.plot(self.r, self.density.mean(axis=0), linewidth=2)

        plt.ylabel('Density (count / nm$^3$)', fontsize=14)
        plt.xlabel('Distance from pore center (nm)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.title('%s' % self.residue.name)
        plt.tight_layout()

        plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    sys = System(args.gro, args.traj, args.residue, args.build_monomer_residue, begin=args.begin, end=args.end,
                 skip=args.skip)

    sys.radial_distribution_function(bins=args.bins)
    sys.bootstrap(args.nboot)
    sys.plot()
