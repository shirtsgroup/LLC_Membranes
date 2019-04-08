#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from LLC_Membranes.analysis import coordination_number, hbonds
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate coordination number')

    # Trajectory Control
    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file. Make sure '
                        'molecules are whole')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Start frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Last frame')
    parser.add_argument('-skip', default=1, type=int, help='Include every skip frames in calculation')

    parser.add_argument('-type', '--type', default='coord', help='Type of lifetime to calculate. Choose `hbonds` to'
                                                                 'calculate hydrogen bond lifetimes. Chosee `coord` to'
                                                                 'calculate lifetime of coordinated molecules')

    # Choose residues

    parser.add_argument('-r', '--residues', nargs='+', default=None, help='For hbonds: residues to include in h-bond '
                        'search. For coordination: Residue to calculate coordination number with respect to.')
    parser.add_argument('-a', '--atoms', action='append', nargs='+', help='Atoms to include from residue. For hbonds, '
                        'each list of atoms must be passed with a separate -a flag for each residue in args.residues')

    # Coordination-specific arguments

    parser.add_argument('-rc', '--coordinated_residue', default=None, help='Name of residue coordinated to residue_')
    parser.add_argument('-ac', '--coordinated_atoms', default=None, nargs='+', help='Name of residue coordinate to '
                        'residue.')
    parser.add_argument('-ta', '--atype', default=None, help='Element name of atoms of which you want coordination '
                                                             'number.')
    parser.add_argument('-tc', '--coordinated_type', default=None, help='Element name of coordinated atoms')
    parser.add_argument('-cut', default=0.25, type=float, help='Maximum distance between pairs where they are '
                                                               'considered coordinated (nm).')

    parser.add_argument('-bins', nargs='+', default=100, type=int, help='Integer or array of bin values. If more than'
                        'one value is used, order the inputs according to the order given in args.axis')
    parser.add_argument('-nocom', action='store_true', help='Calculate coordination based on pairwise distance between'
                                                            'all atoms, even those that are a part of the same residue.')

    # Hbond-specific arguments

    parser.add_argument('-d', '--distance', default=.35, type=float, help='Maximum distance between acceptor and donor'
                        ' atoms. Same as `cut` but separated for now to keep defaults different.')
    parser.add_argument('-angle', '--angle_cut', default=30, type=float,
                        help='Maximum DHA angle to be considered an H-bond')
    parser.add_argument('-acc', '--acceptors', default=False, nargs='+', help='If you only want hbonds with acceptor '
                        'atoms of a residue, use this flag. 1 for this condition, otherwise 0. Input as a list '
                        'of 0 and 1 in the same order as args.residues')
    parser.add_argument('-donors', '--donors', default=False, nargs='+', help='If you only want hbonds with donor atoms'
                        ' of a residue, use this flag. 1 for this condition, otherwise 0. Input as a list of 0 and 1 in'
                        ' the same order as args.residues')

    # saving / loading
    parser.add_argument('-s', '--savename', default='lifetime.pl', help='Name under which to save object.')
    parser.add_argument('-l', '--load', action="store_true", help='Load previously saved object.')

    # display options
    parser.add_argument('-noshow', '--noshow', action="store_true", help='Do not display plots at end')

    return parser


class CoordinationLifetime(coordination_number.System):

    def __init__(self, traj, gro, residue=None, coordinated_residue=None, atoms=None, coordinated_atoms=None, type=None,
                 ctype=None, begin=0, end=-1, skip=1, com=True, cut=0.25):
        """ Narrow system down to groups of interest

        :param traj: Name of GROMACS trajectory (.xtc or .trr)
        :param gro: Name of GROMACS coordinate file (.gro)
        :param residue: Name of residue to include in calculation
        :param coordinated_residue: Name of residue whose coordination we are interested in
        :param atoms: Specify names of atoms to include. All else will be excluded
        :param coordinated_atoms: Specify names of coordinated atoms to include. All else will be excluded
        :param type: Atom types to include
        :param ctype: Coordinated atom types to include
        :param begin: first frame to analyze
        :param end: last frame to analyze
        :param skip: number of frames to skip between analysis steps
        :param com: Calculate coordination based on the center of mass position of the selected atom group
        :param cut: cutoff distance for two species to be considered coordinated.
        """

        super().__init__(traj, gro, residue=residue, coordinated_residue=coordinated_residue, atoms=atoms, com=com,
                         coordinated_atoms=coordinated_atoms, type=type, ctype=ctype, begin=begin, end=end, skip=skip)

        self.distance_search(cut=cut)

        # timeseries for main residue or group of atoms (not coordinated set)
        # for example, if residue = NA and there are 400 NA atoms. The following will be of shape (nframes, 400)
        self.coordination_timeseries = np.zeros([self.t.n_frames, self.distances[0].shape[0]])

        self.create_timeseries()

        self.dwell_times = []
        self.dt = (self.t.time[1] - self.t.time[0]) / 1000  # dt in nanoseconds

    def create_timeseries(self):
        """ Process coordination_number output in order to generate timeseries for each solute
        """

        for t in range(self.t.n_frames):
            self.coordination_timeseries[t, self.distances[t].nonzero()[0]] = self.distances[t].nonzero()[1]

    def calculate_lifetimes(self):

        for x in self.coordination_timeseries.T:
            frame = 0
            while frame < x.size:
                if x[frame] != 0:
                    count = 1
                    # while (frame + count) < x.size and x[frame + count] == x[frame]:
                    #     count += 1
                    # allow one 0 in between if the atoms on both sides are the same
                    while (frame + count) < x.size:
                        if x[frame + count] == x[frame]:
                            count += 1
                        elif (frame + count + 1) < x.size and x[frame + count] == 0 and x[frame + count + 1] == x[frame]:
                            count += 1
                        else:
                            break

                    frame += count
                    if frame < x.size:  # don't count the last dwell time since it was not necessarily finished
                        self.dwell_times.append(count * self.dt)
                else:
                    frame += 1

        print("Mean dwell time: %.2f ns" % np.mean(self.dwell_times))
        print("Median Hbond lifetime %.2f ns" % np.median(self.dwell_times))
        print("Maximum dwell time : %.2f ns" % np.max(self.dwell_times))

    def plot_dwell_time_distribution(self, bins=25):
        """ Plot a histogram of the dwell times

        :param bins: number of bins in histogram
        """

        plt.hist(self.dwell_times, bins=bins)
        plt.ylabel("Frequency")
        plt.xlabel("Dwell Time (ns)")
        plt.show()


class HydrogenBondLifetime(hbonds.System):

    def __init__(self, traj, gro, begin=0, end=-1, skip=1):

        super().__init__(traj, gro, begin=begin, end=end, skip=skip)

        self.dt /= 1000  # convert to nanoseconds
        self.unique_hbonds = {}
        self.dwell_times = []
        self.hbond_timeseries = None

    def identify_unique_hbonds(self):
        """
        """

        ndx = 0

        for frame in self.hbonds:
            for hbond in frame.T:
                lst = str([int(x) for x in hbond[:3]])  # make key out of donor, hydrogen and acceptor trio indices
                if lst not in self.unique_hbonds.keys():
                    self.unique_hbonds[lst] = ndx
                    ndx += 1

    def create_timeseries(self):
        """
        """

        self.hbond_timeseries = np.zeros([self.t.n_frames, len(self.unique_hbonds.keys())], dtype=bool)

        for t, frame in enumerate(self.hbonds):
            for hbond in frame.T:
                lst = str([int(x) for x in hbond[:3]])
                self.hbond_timeseries[t, self.unique_hbonds[lst]] = True

    def calculate_lifetimes(self):

        for x in self.hbond_timeseries.T:
            frame = 0
            while frame < x.size:
                if x[frame]:
                    count = 1
                    while (frame + count) < x.size:
                        if x[frame + count]:
                            count += 1
                        elif (frame + count + 1) < x.size and not x[frame + count] and x[frame + count + 1]:
                            count += 1
                        else:
                            break

                    frame += count
                    if frame < x.size:  # don't count the last dwell time since it was not necessarily finished
                        self.dwell_times.append(count * self.dt)
                else:
                    frame += 1

        print("Mean Hbond lifetime: %.2f ns" % np.mean(self.dwell_times))
        print("Median Hbond lifetime %.2f ns" % np.median(self.dwell_times))
        print("Maximum Hbond lifetime : %.2f ns" % np.max(self.dwell_times))

    def plot_dwell_time_distribution(self, bins=25):
        """ Plot a histogram of the dwell times

        :param bins: number of bins in histogram
        """

        plt.hist(self.dwell_times, bins=bins)
        plt.ylabel("Frequency")
        plt.xlabel("Dwell Time (ns)")
        plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.type == 'coord':

        lifetime = CoordinationLifetime(args.traj, args.gro, atoms=args.atoms[0], coordinated_atoms=args.coordinated_atoms,
                                        residue=args.residues[0], coordinated_residue=args.coordinated_residue,
                                        type=args.atype, ctype=args.coordinated_type, begin=args.begin, end=args.end,
                                        skip=args.skip, com=(not args.nocom))

    elif args.type == 'hbonds':

        # Preprocess the arguments
        if not args.acceptors:
            args.acceptors = [False for i in args.residues]
        else:
            args.acceptors = [bool(int(i)) for i in args.acceptors]

        if not args.donors:
            args.donors = [False for i in args.residues]
        else:
            args.donors = [bool(int(i)) for i in args.donors]

        if not args.atoms:
            args.atoms = [['all'] for r in args.residues]  # a default value

        while len(args.atoms) != len(args.residues):
            args.atoms.append(['all'])

        lifetime = HydrogenBondLifetime(args.traj, args.gro, begin=args.begin, end=args.end, skip=args.skip)

        for i, r in enumerate(args.residues):
            lifetime.set_eligible(r, args.atoms[i], acceptor_only=args.acceptors[i], donors_only=args.donors[i])

        lifetime.identify_hbonds(args.distance, args.angle_cut)
        lifetime.identify_unique_hbonds()
        lifetime.create_timeseries()

    else:

        import sys
        sys.exit('Please choose a type of lifetime to measure (hbonds or coord)')

    lifetime.calculate_lifetimes()

    if not args.noshow:
        lifetime.plot_dwell_time_distribution()
