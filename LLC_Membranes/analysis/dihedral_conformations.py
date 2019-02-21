#!/usr/bin/env python

import argparse
import os
import numpy as np
import mdtraj as md
from LLC_Membranes.analysis import detect_peaks
from LLC_Membranes.llclib import file_rw, physical, topology
import matplotlib.pyplot as plt
import itertools
import sys

"""
Calculate the probabilities for all combinations of chosen dihedral angles
"""


def initialize():

    parser = argparse.ArgumentParser(description='Cacluate probabilities of dihedral combinations')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file. Molecules should be'
                                                                                'made whole (gmx trjconv -pbc whole)')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-d', '--dihedrals', nargs='+', action='append', help='Names of atoms, in'
                        'order of their connectivity, that define the dihedral/torsion we are interested in. NOTE: all'
                        'dihedrals of that type will be analyzed')
    parser.add_argument('-r', '--residue', default='HII', type=str, help='Name of residue that dihedral is apart of, as'
                                                                         'it appears in the .gro file')
    parser.add_argument('-exclude', nargs='+', help='Names of atoms to exclude. If that atom appears in a dihedral,'
                                                    'then that dihedral will not be calculated')
    parser.add_argument('-l', '--load', default=False, type=str, help='Load pickled data')
    parser.add_argument('-b', '--bins', default=100, type=int, help='Number of bins')
    parser.add_argument('-p', '--partition', default=False, type=float, help='(Optional) Radial distance from pore '
                        'center used to divide to regions.')
    parser.add_argument('-bm', '--build_monomer', default='NAcarb11V', type=str, help='Name of monomer used to build'
                        'system. Only necessary if using the partition option.')

    return parser


class Dihedrals(object):

    def __init__(self, gro, traj, dihedrals, resname):

        print('Loading trajectory...', flush=True, end='')
        self.t = md.load(traj, top=gro)
        print('Done!')

        self.ndihedrals = len(dihedrals)
        self.resname = resname

        # get distribution of dihedrals
        for i, d in enumerate(dihedrals):

            # get indices of dihedral atoms
            names = []
            ndx = []
            for a in self.t.topology.atoms:
                if a.name in d and a.residue.name == self.resname:
                    names.append(a.name)
                    ndx.append(a.index)

            if i == 0:
                self.nres = len(ndx) // 4
                self.dihedrals = np.zeros([len(dihedrals), self.t.n_frames, self.nres])

            ordered_ndx = np.array([names[:4].index(x) for x in d])

            ndx = np.reshape(ndx, (len(ndx) // 4, 4))[:, ordered_ndx]

            self.dihedrals[i, ...] = md.compute_dihedrals(self.t, ndx) * (180 / np.pi)

            self.dihedrals[i, ...] = np.where(self.dihedrals[i, ...] < 0, self.dihedrals[i, ...] + 360,
                                              self.dihedrals[i, ...])

            # plt.hist(self.dihedrals[i, ...].flatten(), bins=100)
            # plt.show()
            # exit()

        self.bin_edges = []
        self.translate_bins = []  # will hold dictionaries that tranlate bin numbers to their actual values
        self.all_bin_edges = None
        self.dihedral_bins = None  # to hold bin value of each dihedral on each solute
        self.permutations = None
        self.probs = None

        # for calculating radial dependence of conformations
        self.com = None
        self.partition = None
        self.spline = None

    def partition_solutes(self, r, build_monomer, spts=10):
        """ Calculate solute centers of mass, construct spline through pore and then calculate radial distance of
        solute centers of mass from pore center.

        :param r: distance from pore center defining partition
        :param build_monomer: name of monomer used to build system. Will be used to define pore centers.
        :param spts: number of points used to construct spline

        :type r: float
        :type build_monomer: str
        :type spts: int

        """

        residue = topology.Residue(self.resname)
        residue_indices = np.array([a.index for a in self.t.topology.atoms if a.residue.name == self.resname])

        if residue_indices.size == 0:
            sys.exit("No residue %s found" % self.resname)

        residue_atom_names = [a.name for a in self.t.topology.atoms if a.residue.name == self.resname]
        masses = [residue.mass[x] for x in residue_atom_names[:residue.natoms]]

        print('Calculating centers of mass...', end='', flush=True)
        self.com = physical.center_of_mass(self.t.xyz[:, residue_indices, :], masses)
        print('Done!')

        monomer = topology.LC('%s.gro' % build_monomer)
        pore_defining_atoms = monomer.pore_defining_atoms

        pore_atoms = [a.index for a in self.t.topology.atoms if a.name in pore_defining_atoms and a.residue.name in
                      monomer.residues]

        self.spline = physical.trace_pores(self.t.xyz[:, pore_atoms, :], self.t.unitcell_vectors, spts, npores=4,
                                           progress=True)[0]

        self.partition = physical.partition(self.com, self.spline, r, unitcell=self.t.unitcell_vectors, npores=4,
                                            spline=True)

    def build_spline(self, rep='K', frame=-1):
        """ Build the spline into the last frame of the trajectory

        :param pore_centers: output of physical.avg_pore_loc() with spline=True
        :param rep: name of atom to use to represent spline

        :return: 'spline.gro'
        """
        pos = self.t.xyz[frame, ...]

        for i in range(4):
            pos = np.concatenate((pos, self.spline[frame, i, ...]))
        ids = [a.name for a in self.t.topology.atoms]
        res = [a.residue.name for a in self.t.topology.atoms]
        ids += [rep]*self.spline.shape[2]*self.spline.shape[1]
        res += [rep]*self.spline.shape[2]*self.spline.shape[1]

        file_rw.write_gro_pos(pos, 'spline.gro', ucell=self.t.unitcell_vectors[frame, ...], ids=ids, res=res)

    def define_bins(self, nbins=100):
        """ Define the peaks of the dihedral distribution by selecting cross-over points

        :param nbins: Number of cross-over points that will divide dihedral distributions into nbins - 1 sections

        :type nbins: int

        """

        self.all_bin_edges = np.linspace(0, 360, nbins)

        for i in range(self.ndihedrals):

            bins, edges = np.histogram(self.dihedrals[i, ...].flatten(), bins=np.linspace(0, 360, nbins + 1))

            bin_width = edges[1] - edges[0]
            bin_centers = np.array([i + bin_width / 2 for i in edges[:-1]])
            valleys = detect_peaks.detect_peaks(-bins, mpd=25)

            binned = False
            while not binned:

                plt.plot(bin_centers, bins)
                for i in valleys:
                    plt.plot([bin_centers[i], bin_centers[i]], [0, 1.05*np.amax(bins)], '--', color='black')
                plt.ion()
                plt.show()

                user_input = int(input("Enter 1 if you are happy with the bin edge locations: "))
                if user_input != 1:
                    valleys = [int((float(i) / 360)*nbins) for i in input("Enter list of bin edge locations: ").split()]
                    for k, j in enumerate(valleys):
                        if j == nbins:
                            valleys[k] -= 1  # common mistake
                        elif j > nbins or j < 0:
                            sys.exit('Please enter bin edge locations between 0 and 360')
                    plt.clf()
                else:
                    binned = True

            bin_translation = {}
            for i, v in enumerate(valleys[:-1]):
                bin_translation[i] = bin_centers[np.argmax(bins[v:valleys[i + 1]]) + v]

            plt.close()
            self.translate_bins.append(bin_translation)
            self.bin_edges.append(valleys)

    def bin_dihedrals(self):
        """ For each solute, bin each dihedral into the regions defined in define_bins

        :return:
        """

        # can probably flatten everything to speed the following up. but it's already quick enough
        self.dihedral_bins = np.zeros_like(self.dihedrals, dtype=int)

        for r in range(self.nres):
            for d in range(self.ndihedrals):
                self.dihedral_bins[d, :, r] = np.digitize(self.dihedrals[d, :, r],
                                                          self.all_bin_edges[self.bin_edges[d]]) - 1

    def conformational_probabilities(self):
        """ Calculate the number of times the residue exists with a given combination of dihedrals

        :return:
        """

        self.dihedral_bins = self.dihedral_bins.reshape(self.ndihedrals, self.t.n_frames * self.nres)

        bins = []
        for b in self.bin_edges:
            bins.append([i for i in range(len(b) - 1)])

        self.permutations = list(itertools.product(*bins))

        if self.partition is not None:
            self.partition = self.partition.flatten()
            self.probs = np.zeros([2, len(self.permutations)])
            for i, p in enumerate(self.permutations):
                    for db in self.dihedral_bins[:, self.partition].T:
                        if np.array_equal(db, p):
                            self.probs[0, i] += 1
                    for db in self.dihedral_bins[:, ~self.partition].T:
                        if np.array_equal(db, p):
                            self.probs[1, i] += 1

            self.probs /= np.sum(self.probs, axis=1)[:, np.newaxis]

        else:
            self.probs = np.zeros([len(self.permutations)])

            for i, p in enumerate(self.permutations):
                for db in self.dihedral_bins.T:
                    if np.array_equal(db, p):
                        self.probs[i] += 1

            self.probs /= np.sum(self.probs)

    def print_probabilities(self):
        """ Display the probability of each combination of dihedrals from highest to lowest

        Angle printed are measured counter-clockwise

        :return:
        """

        if self.spline is None:
            self.probs = self.probs[np.newaxis, :]

        for p in self.probs:

            headers = ''
            for i in range(self.ndihedrals):
                headers += '{:^7s}'.format('d%s' % (i + 1))

            headers += ' Probability'
            print(headers)

            ordered = np.argsort(p)[::-1]  # list in descending order

            for i in ordered:
                angles = ['{:^7.1f}'.format(self.translate_bins[k][x]) for k, x in enumerate(self.permutations[i])]
                data = ''
                for k in angles:
                    data += k

                data += '{:^12.2f}'.format(p[i])
                print(data)


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.load:
        d = file_rw.load_object(args.load)
    else:
        d = Dihedrals(args.gro, args.traj, args.dihedrals, args.residue)

        if args.partition:
            d.partition_solutes(args.partition, args.build_monomer)
            d.build_spline()

        d.define_bins(nbins=args.bins)
        d.bin_dihedrals()
        d.conformational_probabilities()
        file_rw.save_object(d, 'dihedral_conformations.pl')

    d.print_probabilities()
