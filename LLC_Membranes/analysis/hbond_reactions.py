#!/usr/bin/env python

from LLC_Membranes.analysis import hbonds
from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt
import argparse


def initialize():

        parser = argparse.ArgumentParser(description='Detect Hydrogen Bonds')

        parser.add_argument('-t', '--trajectory', type=str, default='PR_whole.xtc', help='Name of Gromacs trajectory '
                                                                                          'file (.xtc or .trr)')
        parser.add_argument('-g', '--gro', type=str, default='em.gro', help='Name of Gromacs coordinate file (.gro)')
        parser.add_argument('-d', '--distance', type=float, default=0.3)
        parser.add_argument('-a', '--angle', type=float, default=30)
        parser.add_argument('-r', '--residues', default=('BOC', 'LAL'), help='Names of residues that are involved in'
                                                                             'hydrogen bonding.')
        parser.add_argument('-atoms', '--atoms', action='append', nargs='+', help='Names of atoms involved in hbonding '
                                                                                  'for each residue')
        parser.add_argument('-l', '--load', action="store_true", help='Load pickled data file if it exists')
        parser.add_argument('-show', '--show', action="store_true", help='Show plot at end.')

        return parser


class HbondReactions(hbonds.System):

    def __init__(self, traj, gro):

        super().__init__(traj, gro)

        self.terminated_donors = []
        self.terminated_acceptors = []

    def count_reactions(self, donors=('O1'), acceptors=('N1')):
        """ Count reactions cumulatively over time. Once an hbond forms, it is assumed to have undergone an acid-base
        reaction and the two head groups can no longer react.

        :return:
        """

        for t in range(self.t.n_frames):

            self.terminated_donors.append([])
            self.terminated_acceptors.append([])
            if t > 0:
                self.terminated_donors[t] += self.terminated_donors[t - 1]
                self.terminated_acceptors[t] += self.terminated_acceptors[t - 1]

            for hbond in self.hbonds[t].T:
                a, d = hbond[2], hbond[0]
                if self.names[int(a)] in acceptors and self.names[int(d)] in donors:
                    if a not in self.terminated_acceptors[t] and d not in self.terminated_donors[t]:
                        self.terminated_acceptors[t].append(a)
                        self.terminated_donors[t].append(d)
        #
        # for t in range(self.t.n_frames):
        #
        #     self.terminated_donors.append([])
        #     self.terminated_acceptors.append([])
        #
        #     frame = t
        #     if t == 0:
        #         frame = 1
        #     else:
        #         self.terminated_donors[frame] += self.terminated_donors[frame - 1]
        #
        #     for d, a in self.hbonds[t][[0, 2], :].T:
        #
        #         if d not in self.terminated_donors[frame - 1] and a not in self.terminated_acceptors[frame - 1]:
        #
        #             if self.names[int(d)] in donors and self.names[int(a)] in acceptors:
        #
        #                 self.terminated_donors[t].append(d)
        #                 self.terminated_acceptors[t].append(a)
        #
        #     print(self.terminated_donors[t])
        #     if t == 1:
        #         exit()

    def plot_reaction(self, save=True, show=False, label=None):

        nacceptors = len([x for x in self.D if self.names[x] == 'N1'])
        ndonors = len([x for x in self.D if self.names[x] == 'O1'])

        total_possible = min(ndonors, nacceptors)

        # I think something is wrong if using np.unique is necessary
        # nreacted = [len(np.unique(x)) for x in self.terminated_donors]  # terminated donors and acceptors should be same length

        nreacted = [len(x) for x in self.terminated_donors]
        fraction_reacted = [100 * x / total_possible for x in nreacted]

        plt.plot(self.t.time / 1000, fraction_reacted, lw=2, label=label)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Percent reacted', fontsize=14)
        plt.tight_layout()

        if save:
            plt.savefig('hbonds.pdf')

        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.load:
            print('Loading pickled object...', end='')
            sys = file_rw.load_object('hbonds_d%.2f_a%d.pl' % (10*args.distance, int(args.angle)))
            print('Done')
    else:
            sys = HbondReactions(args.trajectory, args.gro)
            residues = args.residues
            # atoms = {'LAL': ['O1', 'H3'], 'BOC': ['N1']}

            for i, r in enumerate(residues):
                    sys.set_eligible(r, args.atoms[i])

            sys.identify_hbonds(args.distance, args.angle)

            file_rw.save_object(sys, 'hbonds_d%.2f_a%d.pl' % (10*args.distance, int(args.angle)))

    sys.count_reactions()
    sys.plot_reaction(show=args.show)
