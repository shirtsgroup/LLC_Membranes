#!/usr/bin/env python

import argparse
import mdtraj as md
from LLC_Membranes.llclib import physical, topology, file_rw
import numpy as np
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate radial distribution of a monomer component')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='traj_whole.xtc', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-r', '--residue', action='append', nargs='+', help='Name of residue whose '
                        'radial distribution function we want to calculate')
    parser.add_argument('-b', '--build_monomer_residue', default='HII', help='Name of monomer residue used to build'
                                                                             'HII phase')
    parser.add_argument('-bins', '--bins', default=50, type=int, help='Number of bins for histogram of distances')
    parser.add_argument('--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('--end', default=-1, type=int, help='Frame to stop doing calculations')
    parser.add_argument('-skip', '--skip', default=1, type=int, help='Analyze every skip frames')
    parser.add_argument('-l', '--load', help='Name of compressed .npz to load')
    parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials')
    parser.add_argument('-spline', '--spline', action="store_true", help='Trace pore centers using a spline')
    parser.add_argument('-spts', '--spline_pts', default=20, type=int, help='Number of points making up the spline of'
                                                                            'each pore')
    parser.add_argument('-atoms', '--atoms', action='append', nargs='+', help='Names of atoms which '
                        'are a part of residue whose radial distribution function we want to calculate.')
    parser.add_argument('-normalize', action="store_true", help='Make the maximum value of the RDFs equal to 1 for '
                                                                'easier visual comparison (and a single y-axis)')
    parser.add_argument('-cut', default=1.5, type=float, help='Largest distance from pore center to include in '
                                                              'calculation')

    return parser


class System(object):

    def __init__(self, gro, traj, residue, monomer, begin=0, end=-1, skip=1, npores=4, atoms=None):
        """ Calculate the radial distribution of residue in a hexagonal phase LLC Membrane

        :param gro: Coordinate file (.gro or .pdb)
        :param traj: Trajectory file (.trr or .xtc)
        :param residue: Name of residue whose rdf will be calculated
        :param monomer: Name of monomer used to build LLC Membrane
        :param begin: First frame index
        :param end: Last frame index
        :param skip: Skip every 'skip' frames
        :param npores: Number of pores
        :param atoms: Calculate RDF of atoms specified here which are a part of residue

        :type gro: str
        :type traj: str
        :type residue: str
        :type monomer: str
        :type begin: int
        :type end: int
        :type skip: int
        :type npores: int
        :type atoms: list
        """

        self.t = md.load(traj, top=gro)[begin:end:skip]
        self.box = self.t.unitcell_vectors
        self.npores = npores

        if residue == 'SOL':  # workaround for mdtraj
            residue = 'HOH'

        self.residue = topology.Residue(residue)
        self.monomer = topology.LC('%s.gro' % monomer)

        if atoms is not None and 'all' not in atoms:
            res = [a.index for a in self.t.topology.atoms if a.residue.name == residue and a.name in atoms]
            mass = [self.residue.mass[v] for v in self.residue.mass.keys() if v in atoms]
        else:
            res = [a.index for a in self.t.topology.atoms if a.residue.name == residue]
            mass = [v for v in self.residue.mass.values()]

        self.com = physical.center_of_mass(self.t.xyz[:, res, :], mass)

        self.r = None
        self.density = None
        self.bootstraps = None
        self.errorbars = None

    def build_com(self, rep='K'):
        """ Build system with COM of residue plotted in order to confirm that this script is working

        :param rep: name of atom to use to represent spline

        :type rep: str

        :return: structure file, 'spline.gro'
        """

        pos = np.concatenate((self.t.xyz[-1, ...], self.com[-1, ...]))

        ids = [a.name for a in self.t.topology.atoms]
        res = [a.residue.name for a in self.t.topology.atoms]
        ids += [rep] * self.com.shape[1]
        res += [rep] * self.com.shape[1]

        file_rw.write_gro_pos(pos, 'com.gro', ucell=self.box[-1, ...], ids=ids, res=res)

    def radial_distribution_function(self, bins=50, cut=1.5, spline=False, progress=True, npts_spline=10):
        """ Calculate the radial distribution function based on xy distance of solute center of mass from pore center

        :param bins: number of bins in histogram of radial distances
        :param cut: largest distance to include in radial distance histogram
        :param spline: locate pore centers as a function of z. Recommended. Slower, but more accurate
        :param progress: Show progress bar while generating spline
        :param npts_spline: Number of points making up spline tr each pore
        :return:
        """

        self.r = np.zeros([bins])
        self.density = np.zeros([self.t.n_frames, bins])

        pore_defining_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms
                               and a.residue.name == self.monomer.name]

        if spline:
            print('Generating spline through each pore...')

        pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, pore_defining_atoms, :], self.box, spline=spline,
                                             progress=progress, npts=npts_spline)

        if spline:
            self.build_spline(pore_centers)  # to check that the spline was constructed properly
            print('Calculating component density')

        self.r, self.density = physical.compdensity(self.com, pore_centers, self.t.unitcell_vectors,
                                                    nbins=bins, spline=spline, cut=cut)

    def build_spline(self, pore_centers, rep='K'):
        """ Build the spline into the last frame of the trajectory

        :param pore_centers: output of physical.avg_pore_loc() with spline=True
        :param rep: name of atom to use to represent spline

        :return: 'spline.gro'
        """
        pos = self.t.xyz[-1, ...]

        for i in range(4):
            pos = np.concatenate((pos, pore_centers[-1, i, ...]))
        ids = [a.name for a in self.t.topology.atoms]
        res = [a.residue.name for a in self.t.topology.atoms]
        ids += [rep]*pore_centers.shape[2]*pore_centers.shape[1]
        res += [rep]*pore_centers.shape[2]*pore_centers.shape[1]

        file_rw.write_gro_pos(pos, 'spline.gro', ucell=self.box[-1, ...], ids=ids, res=res)

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

    def plot(self, show=False, normalize=False, save=True, savename='rdf.pdf'):

        mean = self.density.mean(axis=0)
        if normalize:
            maximum = np.max(mean)
            mean /= maximum

        if self.bootstraps is not None:
            if normalize:
                self.errorbars /= maximum
            plt.plot(self.r, mean, linewidth=2, label=self.residue.name)
            #np.savez_compressed('NA_rdf.npz', r=self.r, mean=mean)
            plt.fill_between(self.r, self.errorbars[1, :] + mean, mean - self.errorbars[0, :], alpha=0.7)
        else:
            plt.plot(self.r, self.density.mean(axis=0), linewidth=2)

        plt.ylabel('Density (count / nm$^3$)', fontsize=14)
        plt.xlabel('Distance from pore center (nm)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        #plt.title('%s' % self.residue.name)
        plt.legend()
        plt.tight_layout()

        if show:
            plt.show()

        if save:
            plt.savefig(savename)


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.atoms is None:
        args.atoms = [['all']]

    rdfs = []
    for i, r in enumerate(args.residue):

        status = 'Calculating RDF of residue %s' % r[0]
        if args.atoms[i] is not None and args.atoms[i][0] != 'all':
            status += ' restricted to atoms'
            for a in args.atoms[i]:
                status += ' %s' % a
        print(status)

        rdfs.append(System(args.gro, args.traj, r[0], args.build_monomer_residue, begin=args.begin,
                           end=args.end, skip=args.skip, atoms=args.atoms[i]))
        rdfs[i].radial_distribution_function(bins=args.bins, spline=args.spline, npts_spline=args.spline_pts,
                                             cut=args.cut)
        rdfs[i].bootstrap(args.nboot)
        rdfs[i].plot(show=False, normalize=True, save=True)

    #plt.show()