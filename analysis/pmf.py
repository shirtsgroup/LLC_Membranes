#!/usr/bin/env python

import mdtraj as md
import numpy as np
import pymbar
import argparse
from LLC_Membranes.setup.residue_topology import Residue
from LLC_Membranes.analysis import Atom_props
import matplotlib.pyplot as plt
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-f', '--basename', default='umbrella', type=str, help='Name of umbrella simulations where'
                        'i is the state number indexed starting from 0. Assume that all files follow suit. i.e. in the'
                        ' default case, there exists umbrella_i_pullx.xvg and umbrella_i.gro')
    parser.add_argument('-i', '--initial', default='em', type=str, help='Basename corresponding to initial '
                        'configurations used to run umbrella simulations. em_i.gro')
    parser.add_argument('-N', '--nstates', default=12, type=int, help='Number of umbrella simulations')
    parser.add_argument('-r', '--residue', default='ETH', type=str, help='Name of residue whose center of mass was '
                        'constrained')
    parser.add_argument('-c', '--centers', default=None, help='File containing the locations where the center of mass'
                        'was orginally restrained in a column. If None, this script will generate the list for you '
                        'using the initial .gro files')
    parser.add_argument('-k', '--spring_constant', default=1000, type=float, help='spring constant for harmonic'
                                                                                  'restraint (kJ mol^-1 nm^-2)')
    parser.add_argument('-T', '--temperature', default=300, type=int, help='Temperature of simulation (K)')
    parser.add_argument('-b', '--bins', default=50, type=int, help='Number of bins for calculating PMF')

    args = parser.parse_args()

    return args

# Not needed functionality. Delete once everything is working.
# def center_of_mass(positions, res):
#     """
#     :param positions: position of each atom
#     :param res: residue object generated from residue_topology.Residue()
#     :return: center of mass of residue in coordinate file
#     """
#
#     if positions.shape[0] != res.natoms:
#         print('Number of positions does not equal number of atoms in residue. Check your work.')
#         exit()
#
#     com = np.zeros([3])
#     for i in range(res.natoms):
#         com += positions[i, :] * res.masses[i]  # weight each position by the atom's atomic weight
#
#     com /= res.mw
#
#     return com

# if not centers:
#     print('No file with was given with initial COM positions, so I will generate one for you ...',
#           flush=True, end='')
#     t = md.load('%s_0.gro' % self.basename)  # assumes basename_0.gro exists
#     self.res = Residue(residue)
#     self.res_ndx = [a.index for a in t.topology.atoms if a.residue.name == residue]
#     self.nres = len(self.res_ndx) // self.res.natoms
#     self.calculate_initial_com(initial)
#     print('Done. You are welcome. Next time run this script with -c centers.dat unless you want me to do this '
#           'calculation again. Even if this was quick, I prefer that I only have to do this once because I am a '
#           'lazy computer')
# else:
#     with open(centers, 'r') as f:
#         c = []
#         for line in f:
#             c.append(line)
#
#     self.nres = len(c[0].split())
#     self.centers = np.zeros([self.nres, self.K])
#     for i in range(self.K):
#         self.centers[:, i] = c[i].split()

# def calculate_initial_com(self, initial):
#
#     self.centers = np.zeros([self.nres, self.K])  # calculate COM for each residue separately
#     for i in range(self.K):
#         pos = md.load('%s_%d.gro' % (initial, i)).xyz[0, self.res_ndx, :]
#         for j in range(self.nres):
#             self.centers[j, i] = center_of_mass(pos[j*self.res.natoms:(j+1)*self.res.natoms, :], self.res)[2]
#
#     with open('centers.dat', 'w') as f:
#         for i in range(self.K):
#             c = [x for x in self.centers[:, i]]
#             f.write('{:1.3f} {:1.3f} {:1.3f} {:1.3f}\n'.format(c[0], c[1], c[2], c[3]))

# ndx = []
# with open('index.ndx', 'r') as f:
#     for line in f:
#         ndx.append(line)
#
# membrane = 0
# while ndx[membrane].count('[ membrane ]') == 0:
#     membrane += 1
# membrane += 1
#
# indices = []
# while ndx[membrane] != '\n':
#     data = ndx[membrane].split()
#     for i in data:
#         indices.append(int(i))
#     membrane += 1
#
# refcom = np.zeros([12, 1001])
# for i in range(12):
#
#     t = md.load('umbrella_%d.trr' % i, top='umbrella_%d.gro' % i)
#     pos = t.xyz[:, indices, 2]  # z positions of all atoms within reference com index group
#     w = np.array([Atom_props.mass[a.name] for a in t.topology.atoms if a.index in indices])  # weights for com measure
#     refcom[i, :] = np.sum((pos * w), axis=1) / sum(w)
#
# with open('refcom.dat', 'w') as f:
#     for i in range(1001):
#         data = refcom[:, i]
#         for d in data:
#             f.write('%1.3f ' % d)
#         f.write('\n')


class Umbrellas(object):

    def __init__(self, K, basename, residue, k, T, dim='z', centers=None, initial=None):
        """
        :param K: number of umbrella simulations
        :param basename: Name of umbrella simulations where i is the state number indexed starting from 0. Assumes that
                        all files follow suit. i.e. there exists umbrella_i_pullx.xvg and umbrella_i.gro
        :param residue: Name of residue whose center of mass was constrained
        :param k : spring constant for harmonic restraint
        """

        kB = 1.381e-23 * 6.022e23 / 1000.0  # Boltzmann constant in kJ/mol/K
        self.beta = 1 / (kB * T)
        self.spring_constant = k  # (kJ mol^-1 nm^-2)
        self.K = K
        self.centers = np.zeros([self.K])
        self.basename = basename
        self.uncorrelated_samples = []
        self.com_bin = []
        self.histo = []
        self.bin_centers = []  # for plotting
        self.results = None

        self.dimension = []
        for d in 'xyz':
            if d in dim:
                self.dimension.append(d)

        # read pullx files and store com location with each frame
        # do the first one outside of the loop so we can intialize com_dist array to be the right size
        pullx = []
        with open('%s_1_pullx.xvg' % basename, 'r') as f:
            for line in f:
                if line[0] != '#' and line[0] != '@':
                    pullx.append(line)

        # first column is time. Two columns per pull group (absolute distance, z-component distance)
        self.nres = (len(pullx[0].split()) - 1) // 2

        self.com_dist = np.zeros([self.K, self.nres, len(pullx)])  # initialize

        c = []
        with open('centers.txt', 'r') as f:
            for line in f:
                c.append(line)

        initial_com = np.zeros([self.K, self.nres])
        for i in range(1, self.K + 1):
            initial_com[i - 1, :] = c[i].split()

        # ref = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
        # t = md.load('long_1.trr', top='long_1.gro')
        # ref_atoms = [a.index for a in t.topology.atoms if a.name in ref and 1369 <= a.index < 2055]
        # eth = [a.index for a in t.topology.atoms if a.residue.name == 'ETH']
        # eth_names = [a.name for a in t.topology.atoms if a.residue.name == 'ETH']
        # w = [Atom_props.mass[a] for a in eth_names]
        # eth = eth[:9]
        # w = w[:9]
        #
        # ref_com = np.zeros([t.n_frames])
        # eth_com = np.zeros([t.n_frames])
        # for i in range(t.n_frames):
        #     ref_com[i] = np.mean(t.xyz[i, ref_atoms, 2])
        #     for j in range(9):
        #         eth_com[i] += w[j]*t.xyz[i, eth[j], 2]
        #     eth_com[i] /= sum(w)
        #
        # plt.hist(np.abs(eth_com - ref_com), bins=25)
        # plt.show()
        #
        # exit()

        print('Reading pullx files...', end='', flush=True)
        for i in range(self.K):
            if i != 0:
                pullx = []
                with open('%s_%d_pullx.xvg' % (basename, i + 1), 'r') as f:
                    for line in f:
                        if line[0] != '#' and line[0] != '@':
                            pullx.append(line)
            for j in range(self.com_dist.shape[2]):
                self.com_dist[i, :, j] = pullx[j].split()[2::2]  # extract COM dZ

        self.com_dist += initial_com[:, :, None]  # add initial com location to properly space apart histograms

        colors = ['aqua', 'blue', 'coral', 'crimson', 'darkgreen', 'gold', 'lavender', 'magenta', 'orangered', 'plum',
                  'teal', 'violet']
        self.centers = self.com_dist[:, :, 0]

        for i in range(self.K):
            plt.hist(self.com_dist[i, 0, :], bins=50, color=colors[i])
            plt.plot([self.centers[i, 0], self.centers[i, 0]], [0, 10000], '--', color=colors[i])
            # print(np.mean(self.com_dist[i, 0, :]))
            # print(np.std(self.com_dist[i, 0, :]))

        plt.show()
        exit()

        self.N = np.zeros([self.K, self.nres], dtype=int)  # number of uncorrelated frames for each trajectory
        self.u_kn = np.zeros_like(self.com_dist)
        self.u_kln = np.zeros([self.nres, self.K, self.K, len(pullx)])

    def extract_uncorrelated_samples(self, t=None):
        """
        :param t: correlation time. Set the correlation time manually (in terms of frames between uncorrelated samples)
        or let pymbar figure it out with the timeseries module by default (None)
        """

        g = t
        print('Extracting uncorrelated samples...', end='', flush=True)
        for u in range(self.K):
            for r in range(self.nres):
                if not t:
                    g = pymbar.timeseries.statisticalInefficiency(self.com_dist[u, r, :], fft=True)
                indices = pymbar.timeseries.subsampleCorrelatedData(self.com_dist[u, r, :], g=g)
                self.N[u, r] = len(indices)
                self.com_dist[u, r, :len(indices)] = self.com_dist[u, r, indices]
        print('Done!')

    def calculate_PMF(self, nbins):

        #dz = self.com_dist - self.centers[..., np.newaxis]  # deviation from restrained position

        self.com_bin = np.zeros_like(self.com_dist, dtype=int)  # [umbrella, residue, config]
        self.histo = np.zeros([self.nres, self.K, nbins], dtype=int)
        self.bin_centers = np.zeros([self.nres, nbins])
        self.results = []

        for i in range(self.nres):
            print('Calculating PMF for Residue %d' % i, flush=True)
            # left edge of bins
            maxi = np.amax(self.com_dist[:, i, :])
            mini = np.amin(self.com_dist[:, i, :])
            delta = (maxi - mini) / nbins
            bins = np.linspace(mini, maxi - delta, nbins)  # some reasonable bounds
            self.bin_centers[i, :] = bins + 0.5*delta
            bins += delta  # for proper output from np.digitize
            for k in range(self.K):
                for n in range(self.N[k, i]):
                    dz = self.com_dist[k, i, n] - self.centers[:, i]
                    self.u_kln[i, k, :, n] = self.u_kn[k, i, n] + 0.5 * self.beta * self.spring_constant * dz**2
                    self.com_bin[k, i, n] = np.digitize(self.com_dist[k, i, n], bins)
                self.histo[i, k, :] = [np.sum(self.com_bin[k, i, :self.N[k, i]] == a) for a in range(nbins)]
            mbar = pymbar.MBAR(self.u_kln[i, ...], self.N[:, i])
            self.results.append(mbar.computePMF(self.u_kn[:, i, :], self.com_bin[:, i, :], nbins))

    def plot_histograms(self, show=True, save=False, savename='histograms.png'):

        for i in range(self.nres):
            plt.figure(i + 1)
            for k in range(self.K):
                plt.plot(self.bin_centers[i, :], self.histo[i, k, :])

        plt.tight_layout()
        if save:
            plt.savefig(savename)
        if show:
            plt.show()

    def plot_PMF(self, show=True, save=False, savename='pmf.png'):

        n = int(np.ceil(self.nres ** 0.5))
        fig, ax = plt.subplots(n, n)

        for i in range(n):
            for j in range(n):
                conf = i * n + j
                if conf < self.nres:
                    ax[i, j].errorbar(self.bin_centers[conf, :], self.results[conf]['f_i'], yerr=self.results[conf]['df_i'])
                    ax[i, j].set_title('Residue %d' % conf)
                    ax[i, j].set_ylabel('Free Energy (kJ/mol)')
                    ax[i, j].set_xlabel('Distance from COM reference (nm)')

        # plt.tight_layout()
        if save:
            plt.savefig(savename)
        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize()

    # initialize umbrella data
    u = Umbrellas(args.nstates, args.basename, args.residue, args.spring_constant, args.temperature,
                  initial=args.initial, centers=args.centers)
    u.extract_uncorrelated_samples(t = 20)
    u.calculate_PMF(args.bins)
    # u.plot_histograms()
    u.plot_PMF()
