#!/usr/bin/env python

import argparse
import mdtraj as md
from LLC_Membranes.llclib import physical, topology, file_rw
import numpy as np
import matplotlib.pyplot as plt
import sys


def initialize():

    parser = argparse.ArgumentParser(description='Calculate radial distribution function for components of a hexagonal'
                                                 ' phase LLC membrane system from the command line')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='traj_whole.xtc', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-r', '--residue', action='append', nargs='+', help='Name of residue whose '
                        'radial distribution function we want to calculate')
    parser.add_argument('-b', '--build_monomer_residue', default='NAcarb11V', help='Name of monomer residue used to build'
                                                                             'HII phase')
    parser.add_argument('-bins', '--bins', default=50, type=int, help='Number of bins for histogram of distances')
    parser.add_argument('--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('--end', default=-1, type=int, help='Frame to stop doing calculations')
    parser.add_argument('-skip', '--skip', default=1, type=int, help='Analyze every skip frames')
    parser.add_argument('-l', '--load', help='Name of compressed .npz to load')
    parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials')
    parser.add_argument('-spline', '--spline', action="store_true", help='Trace pore centers using a spline')
    parser.add_argument('-spts', '--spline_pts', default=10, type=int, help='Number of points making up the spline of'
                                                                            'each pore')
    parser.add_argument('-atoms', '--atoms', action='append', nargs='+', help='Names of atom(s) which '
                        'is/are a part of residue whose radial distribution function we want to calculate.')
    parser.add_argument('-normalize', action="store_true", help='Make the maximum value of the RDFs equal to 1 for '
                                                                'easier visual comparison (and a single y-axis)')
    parser.add_argument('-cut', default=1.5, type=float, help='Largest distance from pore center to include in '
                                                              'calculation')

    return parser


def grps(name):
    """ Return names of atoms making up specialized groups. This is specifically for NaGA3C11.

    TODO: incorporate these groups as annotations (maybe)

    :param name: specialized group names. Valid options are 'head groups' and 'tails'

    :type name: str

    :return: atom names that constitute specialized groups
    """

    if name == 'head groups':

        return ['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O3', 'O4']

    elif name == 'tails':

        return ['C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21',
                'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35',
                'C36', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48']


class System(object):
    """ Calculate the radial distribution of a residue, atom or group of atoms in a hexagonal phase LLC Membrane
    """

    def __init__(self, gro, traj, residue, monomer, begin=0, end=-1, skip=1, npores=4, atoms=None, com=True):
        """ Restrict the system to desired components. Calculate centers of mass of components if desired

        :param gro: GROMACS coordinate file (.gro or .pdb)
        :param traj: GROMACS trajectory file (.trr or .xtc)
        :param residue: Name of residue to which the desired component belongs. If no other options are specified the \
        radial distribution of this residue's center of mass will be calculated
        :param monomer: Name of monomer used to build LLC Membrane
        :param begin: Index of first frame of trajectory to be used in analysis
        :param end: Index of last frame of trajectory to be used in analysis
        :param skip: Do not include every 'skip' frames in the analysis
        :param npores: Number of pores in the unit cell (Currently only works for 4 pores)
        :param atoms: Restrict residue to atoms named here
        :param com: Calculate density based on center of mass of residue or group of atoms

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

        print('Loading Trajectory...', end='', flush=True)
        self.t = md.load(traj, top=gro)[begin:end:skip]
        print('Done!')
        self.box = self.t.unitcell_vectors
        self.npores = npores

        names = topology.fix_names(gro)  # rename atoms because mdtraj screws it up in some cases.
        for i, a in enumerate(self.t.topology.atoms):
            a.name = names[i]

        if residue == 'SOL':  # workaround for mdtraj
            residue = 'HOH'

        self.residue = topology.Residue(residue)
        self.monomer = topology.LC('%s' % monomer)

        if atoms is not None and 'all' not in atoms:
            res = [a.index for a in self.t.topology.atoms if a.residue.name == residue and a.name in atoms]
            if com:
                mass = [self.residue.mass[v] for v in self.residue.mass.keys() if v in atoms]
        else:
            res = [a.index for a in self.t.topology.atoms if a.residue.name == residue]
            if com:
                mass = [v for v in self.residue.mass.values()]

        if com:
            print('Calculating centers of mass...', flush=True, end='')
            self.com = physical.center_of_mass(self.t.xyz[:, res, :], mass)
            print('Done!')
        else:
            self.com = self.t.xyz[:, res, :]

        self.r = None
        self.density = None
        self.bootstraps = None
        self.errorbars = None
        self.radial_distances = None

    def build_com(self, rep='K'):
        """ Build system with the center of mass of the residue or group of atoms plotted in order to confirm that this script is working

        :param rep: name of atom to use to represent spline

        :type rep: str

        :return: 'com.gro'
        """

        pos = np.concatenate((self.t.xyz[-1, ...], self.com[-1, ...]))

        ids = [a.name for a in self.t.topology.atoms]
        res = [a.residue.name for a in self.t.topology.atoms]
        ids += [rep] * self.com.shape[1]
        res += [rep] * self.com.shape[1]

        file_rw.write_gro_pos(pos, 'com.gro', ucell=self.box[-1, ...], ids=ids, res=res)

    def radial_distribution_function(self, bins=50, cut=1.5, spline=False, progress=True, npts_spline=10, build=True):
        """ Calculate the radial distribution function based on xy distance of solute center of mass from pore center

        :param bins: number of bins in histogram of radial distances
        :param cut: largest distance to include in radial distance histogram
        :param spline: locate pore centers as a function of z. Recommended. Slower, but more accurate
        :param progress: Show progress bar while generating spline
        :param npts_spline: Number of points making up spline through each pore
        :param build: Create a .gro file of the last frame with the spline represented as atoms

        :type bins: int
        :type cut: float
        :type spline: bool
        :type progress: bool
        :type npts_spline: int
        :type build: bool
        """

        self.r = np.zeros([bins])
        self.density = np.zeros([self.t.n_frames, bins])

        pore_defining_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms
                               and a.residue.name in self.monomer.residues]

        if not pore_defining_atoms:
            sys.exit('There are no atoms specified which can be used to define the pore centers.\n'
                     'There are a few reasons this might have happened:\n'
                     '1) Did you specify the correct monomer for this system?\n'
                     '2) If so, have you properly annotated the pore defining atoms?\n'
                     '3) Perhaps there is an annotated .itp in LLC_Membranes/top/topologies, but an unannotated .itp '
                     'file for the same monomer present in this directory.\n')

        pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, pore_defining_atoms, :], self.box, spline=spline,
                                             progress=progress, npts=npts_spline)

        if spline and build:
            self.build_spline(pore_centers)  # to check that the spline was constructed properly
            print('Calculating component density')

        self.r, self.density, self.radial_distances = physical.compdensity(self.com, pore_centers,
                                                                           self.t.unitcell_vectors,
                                                                           nbins=bins, spline=spline, cut=cut,
                                                                           radial_distances=True)

    def build_spline(self, pore_centers, rep='K'):
        """ Build the spline into the last frame of the trajectory

        :param pore_centers: output of physical.avg_pore_loc() with spline=True
        :param rep: name of atom to use to represent spline

        :type pore_centers: np.ndarray
        :type rep: str

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
        """ Generate statistics using the bootstrapping technique

        :param nboot: number of bootstrap trials (i.e. number of time data is resampled)

        :type nboot: int
        """

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

    def plot(self, show=False, normalize=False, save=True, savename='rdf.pdf', label=None):
        """ Plot the radial distribution function

        :param show: Display plot
        :param normalize: Divide by all values of RDF by the maximum value
        :param save: Save the plot under savename
        :param savename: Name and format (determined by file extension) under which to save plot
        :param label: Legend label. If None, residue name will be used.

        :type show: bool
        :type normalize: bool
        :type save: bool
        :type savename: str
        :type label: str
        """

        if label is not None:
            lab = label.title()
        else:
            lab = self.residue.name

        mean = self.density.mean(axis=0)
        if normalize:
            maximum = np.max(mean)
            mean /= maximum

        if self.bootstraps is not None:
            if normalize:
                self.errorbars /= maximum
            plt.plot(self.r, mean, linewidth=2, label=lab)
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

    if args.residue is None:
        sys.exit('Please specify a residue to calculate rdf with respect to')  # need better grammar here

    while len(args.atoms) < len(args.residue):
        args.atoms.append(None)

    special_groups = ['head groups', 'tails']

    rdfs = []
    for i, r in enumerate(args.residue):

        com = True

        status = 'Calculating RDF of residue %s' % r[0]
        savename = 'rdf_%s' % r[0]
        if args.atoms[i] is not None and args.atoms[i][0] != 'all':

            if args.atoms[i][0].lower() in special_groups:

                label = args.atoms[i][0].lower()
                status += ' restricted to %s atoms' % label
                savename += '_%s' % label
                args.atoms[i] = grps(label)
                com = False

            else:

                status += ' restricted to atoms'
                savename += '_'
                for a in args.atoms[i]:
                    status += ' %s' % a
                    savename += '%s' % a

                label = None
        else:
            label = None

        print(status)

        rdfs.append(System(args.gro, args.traj, r[0], args.build_monomer_residue, begin=args.begin,
                           end=args.end, skip=args.skip, atoms=args.atoms[i], com=com))
        rdfs[i].radial_distribution_function(bins=args.bins, spline=args.spline, npts_spline=args.spline_pts,
                                             cut=args.cut)
        rdfs[i].bootstrap(args.nboot)
        rdfs[i].plot(show=False, normalize=args.normalize, save=True, label=label)

        file_rw.save_object(rdfs[i], '%s.pl' % savename)

    plt.show()
