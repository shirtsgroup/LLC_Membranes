#!/usr/bin/env python

import argparse
from LLC_Membranes.analysis import hbonds
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='wiggle.trr', type=str, help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate file')
    parser.add_argument('-p', '--top', default='topol.top', type=str, help='Gromacs topology file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='End frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Start frame')
    parser.add_argument('-x', '--exclude_water', action='store_true', help='Exclude water while searching for hbonds')
    parser.add_argument('-r', '--residues', nargs='+', default=['HII'], help='Residues to include in h-bond search. '
                        'Water is automatically included if you do not specify the -x option')
    parser.add_argument('-a', '--atoms', action='append', nargs='+', help='Atoms to include for each '
                        'residue. Each list of atoms must be passed with a separate -a flag for each residue in'
                        'args.residues')
    parser.add_argument('-d', '--distance', default=.3, help='Maximum distance between acceptor and donor atoms')
    parser.add_argument('-angle', '--angle_cut', default=20, help='Maximum DHA angle to be considered an H-bond')
    parser.add_argument('-l', '--load', default=False, help='Load pickled hbond trajectory')
    parser.add_argument('-s', '--savename', default='hbonds.pl', help='Name of hbond trajectory file to save')
    parser.add_argument('-nb', '--nboot', default=200, type=int, help='Number of bootstrap trials')
    parser.add_argument('-nbins', '--nbins', default=50, type=int, help='Number of bins in histogram')

    return parser


def wait_times(t, A, lam):

    return A*np.exp(-lam * t)


if __name__ == '__main__':

    args = initialize().parse_args()

    if not args.load:

        # workaround for argparse. If default value is set, it is always included in the list with action='append'
        if not args.atoms:
            args.atoms = [['O3', 'O4']]  # a default value

        sys = hbonds.System(args.traj, args.gro, args.top, begin=args.begin, end=args.end, exclude_water=args.exclude_water)

        for i, r in enumerate(args.residues):
            sys.set_eligible(r, args.atoms[i])

        sys.identify_hbonds(args.distance, args.angle_cut)

        sys.hbond_matrix()

        with open(args.savename, "wb") as f:
            pickle.dump(sys, f)

    else:

        with open(args.load, "rb") as f:
            sys = pickle.load(f)

    # each entry in sys.hbonds is a numpy array of shape (4, nhbonds). Each of the four values is used to describe a
    # single hydrogen bond as follows: [Donor index, Hydrogen index, Acceptor index, angle between DH and HA]

    # The donor-acceptor matrix is of shape [nframes, natoms] where natoms is the number of unique atoms involved in
    # hydrogen bonding. Element [i, j] specifies the index (with respect to all atoms in the system) of the acceptor
    # atom recieving an hbond from donor atom j at frame i. Donor atom j can be converted to the atom index with respect
    # to the entire system by using sys.matrix_to_atom_index[j]

    dwell_time = []
    nonzero = 0
    for atom in range(sys.donor_acceptor_matrix.shape[1]):
        # Don't look at atoms that did appear to have not hydrogen bonded since they are actually just acceptors
        x = sys.donor_acceptor_matrix[:, atom]
        if len(np.nonzero(x)[0]) > 0:
            nonzero += 1
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
                        dwell_time.append(count*sys.dt)
                else:
                    frame += 1

    # bootstrap dwell times
    bootstrap_trials = np.zeros([args.nboot, args.nbins])
    lambda_distribution = np.zeros([args.nboot])
    A = 0
    for b in range(args.nboot):
        trial = np.random.choice(dwell_time, size=len(dwell_time), replace=True)
        bootstrap_trials[b, :], bins = np.histogram(trial, bins=args.nbins, range=(np.min(dwell_time), np.max(dwell_time)))
        p0 = [bootstrap_trials[b, 0], len(dwell_time) / (nonzero * sys.t.time[-1])]
        solp_dwell, cov_x = curve_fit(wait_times, bins[:-1], bootstrap_trials[b, :], p0)
        lambda_distribution[b] = solp_dwell[1]
        A += solp_dwell[0]

    A /= args.nboot

    bin_width = bins[1] - bins[0]
    plt.bar(bins[:-1], bootstrap_trials.mean(axis=0), width=bin_width)
    plt.plot(np.linspace(bins[0], bins[-1], 1000), wait_times(np.linspace(bins[0], bins[-1], 1000), A,
            np.mean(lambda_distribution)), '--', color='black', label='$\lambda_{fit}$ = %.4f $\pm$ %.4f' %
            (np.mean(lambda_distribution), np.std(lambda_distribution)))
    plt.legend(fontsize=14)
    plt.gcf().get_axes()[0].tick_params(labelsize=14)
    plt.xlabel('Hydrogen bond lifetime (ps)', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.tight_layout()
    plt.show()
