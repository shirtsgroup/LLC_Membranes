#!/usr/bin/env python

from LLC_Membranes.analysis import hbonds
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import argparse
import numpy as np


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
        parser.add_argument('-tcl', '--tcl', action="store_true", help='Create a tcl file for viewing reacted groups.')

        return parser


class HbondReactions(hbonds.System):

    def __init__(self, traj, gro):

        super().__init__(traj, gro)

        self.terminated_donors = []
        self.terminated_acceptors = []
        self.fraction_reacted = None

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

    def get_fraction_reacted(self):

        nacceptors = len([x for x in self.D if self.names[x] == 'N1'])
        ndonors = len([x for x in self.D if self.names[x] == 'O1'])

        total_possible = min(ndonors, nacceptors)

        # I think something is wrong if using np.unique is necessary
        # nreacted = [len(np.unique(x)) for x in self.terminated_donors]  # terminated donors and acceptors should be same length

        nreacted = [len(x) for x in self.terminated_donors]
        self.fraction_reacted = [100 * x / total_possible for x in nreacted]

    def plot_reaction(self, save=True, show=False, label=None):

        if self.fraction_reacted is None:
            self.get_fraction_reacted()

        plt.plot(self.t.time / 1000, self.fraction_reacted, lw=2, label=label)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Percent reacted', fontsize=14)
        plt.tight_layout()

        if save:
            plt.savefig('hbonds.pdf')

        if show:
            plt.show()

    def _head_group_indices(self):
        """ Extract the indices of the atoms in each head group that will be colored

        :return hg_indices: dictionary with keys that are residue numbers and values that are lists of indices
        :rtype hg_indices: dict
        """

        # this stuff is hard-coded for now but can be generlized using topology.LC and annoations
        natoms = {'LAL': 148, 'BOC': 147}
        hgatoms = {
            'LAL': ['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O', 'N', 'C50', 'C51', 'C7', 'O1', 'O11', 'H3', 'H80',
                    'H81', 'H79', 'H2', 'H82'],
            'BOC': ['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O', 'N', 'H2', 'C7', 'H3', 'H80', 'C8', 'H81', 'H82',
                    'N1', 'H4', 'H83']}

        res_no = 0
        tot = 0
        res_numbers = []
        residue_names = dict()
        previous_residue = None
        for nr, a in enumerate(self.t.topology.atoms):
            residue = a.residue.name
            if nr == 1:
                previous_residue = residue
            if nr > 0 and (nr - tot) % natoms[previous_residue] == 0:
                residue_names[res_no] = residue
                tot = nr
                #tot += natoms[previous_residue]
                res_no += 1
            res_numbers.append(res_no)
            previous_residue = residue

        residue_names[res_no] = residue

        hg_serial = dict()
        for r in range(res_no + 1):
            resname = residue_names[r]
            serial = []
            # Some brute force
            for a in sys.t.topology.atoms:
                if res_numbers[a.index] == r and a.name in hgatoms[resname]:
                    serial.append(a.index)
            hg_serial[r] = serial

        return hg_serial, res_numbers

    def visualize_reactions(self, frame=-1, name='reacted.tcl'):
        """ Create a tcl file that makes it easy to view reacted an unreacted head groups. Unreacted acid and base
        head groups will be colored red and blue respectively. Reacted groups will be colored gray.

        Currently this is specific to the LAL-BOC system

        :param frame: frame number to visualize
        :param name: name of output tcl file

        :type frame: int
        :type name: str
        """

        acid = np.array([a.index for a in self.t.topology.atoms if a.name == 'O1' and a.residue.name == 'LAL'])
        base = np.array([a.index for a in self.t.topology.atoms if a.name == 'N1' and a.residue.name == 'BOC'])

        hg_serial, residue_numbers = self._head_group_indices()

        reacted_acids = [a for a in acid if a in self.terminated_donors[frame]]
        reacted_bases = [b for b in base if b in self.terminated_acceptors[frame]]
        reacted = reacted_acids + reacted_bases

        unreacted_acids = [a for a in acid if a not in self.terminated_donors[frame]]
        unreacted_bases = [b for b in base if b not in self.terminated_acceptors[frame]]

        with open(name, 'w') as f:

            f.write('color Display Background white\n')

            # Acids
            f.write('mol addrep 0\n')
            f.write('mol modselect 0 0 index')
            for a in unreacted_acids:
                for ndx in hg_serial[residue_numbers[a]]:
                    f.write(' %d' % ndx)
            f.write('\n')
            f.write('mol modcolor 0 0 ColorID 1\n')  # make acids red
            f.write('mol modstyle 0 0 Licorice 0.5 12.0 12.0\n')

            # Bases
            f.write('mol addrep 0\n')
            f.write('mol modselect 1 0 index')
            for b in unreacted_bases:
                for ndx in hg_serial[residue_numbers[b]]:
                    f.write(' %d' % ndx)
            f.write('\n')
            f.write('mol modcolor 1 0 ColorID 0\n')  # make bases blue
            f.write('mol modstyle 1 0 Licorice 0.5 12.0 12.0\n')

            # Unreacted
            f.write('mol addrep 0\n')
            f.write('mol modselect 2 0 index')
            for u in reacted:
                for ndx in hg_serial[residue_numbers[u]]:
                    f.write(' %d' % ndx)
            f.write('\n')
            f.write('mol modcolor 2 0 ColorID 8\n')  # make unreacted white
            f.write('mol modstyle 2 0 Licorice 0.5 12.0 12.0\n')


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

    if args.tcl:
        sys.visualize_reactions()

    sys.plot_reaction(show=args.show)
