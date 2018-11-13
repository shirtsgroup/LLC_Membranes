#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import physical
import pickle


def initialize():

    parser = argparse.ArgumentParser(description='Figure out the weight percent of water in the pores and tails')

    # trajectory control
    parser.add_argument('-t', '--traj', default=False, help='Name of GROMACS trajectory file. This file should be'
                        'preprocessed so everything is the box. For example, use gmx trjconv -ur tric -pbc atom. In the'
                        'event that a box vector crosses through a pore, use shift_box.py first to fix that. Specify'
                                                            'False (default) if you are only looking at a single frame')
    parser.add_argument('-g', '--gro', default='PR.gro', help='Name of GROMACS coordinate file')
    parser.add_argument('-r', '--residue', default='SOL', help='Name of residue whose partition we wish to quantify')
    parser.add_argument('-begin', default=0, type=int, help='First frame to read')
    parser.add_argument('-end', default=-1, type=int, help='Last frame to read')
    parser.add_argument('-skip', default=1, type=int, help='Skip every n frames')

    # define system
    parser.add_argument('-p', '--pore_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms that'
                        'will be used to define the pore region')
    parser.add_argument('-ox', '--tail_oxygen', nargs='+', default=['O5', 'O6', 'O7', 'O8', 'O9', 'O10'], help='Oxygen'
                        'atoms that will be used to define the tail region')
    parser.add_argument('-tr', '--tail_radius', default=0.5, type=float, help='Max distance from tail oxygens a water '
                        'molecule can exist in order to be counted as inside the pore')

    parser.add_argument('-pr', '--pore_radius', default=0.5, type=float, help='Max distance from pore center a water '
                        'molecule can exist in order to be counted as inside the pore')
    parser.add_argument('-b', '--bounds', default=5, type=float, help='Distance from z-center up until which all atoms '
                                                                      'will be included in calculation (nm)')
    parser.add_argument('-natoms', default=137, type=int, help='Number of atoms in monomer residue (not including ions '
                                                               ' if they are separate residues!')

    # save/load options
    parser.add_argument('--savename', default='pore_spline.pl', help='Name of file in which to save System object')
    parser.add_argument('--load', action="store_true")

    parser.add_argument('-boot', '--nboot', default=200, type=int, help='Number of bootstrap trials')
    parser.add_argument('--single_frame', action='store_true', help='Specify this flag in order to analyze a single'
                                                                    '.gro file. No statistics will be generated')

    return parser


class System(object):

    def __init__(self, gro, pore_atoms, residue, traj=False, begin=0, end=-1, skip=1, npores=4):
        """ Define the system and boundaries for pore and tail region

        :param gro: coordinate file
        :param pore_atoms: atoms used to define the pore locations
        :param traj: trajectory file
        :param begin: first frame to include
        :param end: last frame to include
        :param skip: skip every n frames
        :param npores: number of pores. Assumes that atoms are number sequentially by pore
        """

        print('Loading trajectory...', flush=True, end='')
        if traj:
            self.t = md.load(traj, top=args.gro)[begin:end:skip]
        else:
            self.t = md.load(gro)
        print('Done')

        # coordinates and unit cell dimensions
        self.pos = self.t.xyz
        box = self.t.unitcell_vectors
        self.box = [box[0, 0, 0], box[0, 1, 1], box[0, 2, 2], box[0, 0, 1], box[0, 2, 0],
                   box[0, 1, 0], box[0, 0, 2], box[0, 1, 2], box[0, 2, 0]]  # gromacs format
        self.res = np.array([a.residue.name for a in self.t.topology.atoms])  # all of the residues
        self.ids = np.array([a.name for a in self.t.topology.atoms])  # all of the atom names

        # find pore centers
        print('Creating pore splines')
        pore_atoms = [a.index for a in self.t.topology.atoms if a.name in pore_atoms]
        self.pore_spline, self.bin_centers = physical.trace_pores(self.pos[:, pore_atoms, :], self.t.unitcell_vectors,
                                                                  20)

    def plot(self, frame):

        fig, ax = plt.subplots(2, 2, figsize=(10, 10))

        for i in range(4):

            ax1 = ax[i // 2, i % 2]
            ax2 = ax1.twinx()
            spline = self.pore_spline[frame, i, ...]
            bins = self.bin_centers[frame, i, :]

            xrange = (np.amax(spline[:, 0]) - np.amin(spline[:, 0])) / 2
            yrange = (np.amax(spline[:, 1]) - np.amin(spline[:, 1])) / 2

            ax1.plot(bins, spline[:, 0], color='xkcd:blue', linewidth=2)
            ax2.plot(bins, spline[:, 1], color='xkcd:orange', linewidth=2)

            if i % 2 == 0:
                ax1.set_ylabel('$x$-coordinate', fontsize=14, color='xkcd:blue')
            if i % 2 == 1:
                ax2.set_ylabel('$y$-coordinate', fontsize=14, color='xkcd:orange')
            if i // 2 == 1:
                ax1.set_xlabel('$z$-coordinate', fontsize=14)

            # set limits -- give a little white space above and below
            ax1.set_ylim(spline[:, 0].mean() - xrange*2, spline[:, 0].mean() + xrange*2)
            ax2.set_ylim(spline[:, 1].mean() - yrange * 2, spline[:, 1].mean() + yrange * 2)

            # format tick size
            plt.gcf().get_axes()[i].tick_params(labelsize=14)
            ax2.yaxis.set_tick_params(labelsize=14)

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if not args.load:
        sys = System(args.gro, args.pore_atoms, args.residue, traj=args.traj, begin=args.begin, end=args.end,
                     skip=args.skip)

        with open(args.savename, "wb") as f:
            pickle.dump(sys, f)
    else:

        with open(args.savename, "rb") as f:
            sys = pickle.load(f)

    sys.plot(-1)
