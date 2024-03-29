#!/usr/bin/env python

import argparse
import yaml
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from mpl_toolkits.axes_grid1 import AxesGrid
from LLC_Membranes.analysis import hbonds, coordination_number
from LLC_Membranes.llclib import file_rw, physical, topology, timeseries, fitting_functions, sampling
from LLC_Membranes.timeseries.fractional_levy_motion import FLM
from LLC_Membranes.timeseries import flm_sim_params
import sys
import tqdm
from scipy.stats import cauchy, laplace, t, levy_stable
import levy
import itertools
from multiprocessing import Pool
np.set_printoptions(precision=4, suppress=True)
import time as timer
import fbm


def initialize():

    parser = argparse.ArgumentParser(description='Model a continuous time random walk')

    # Pass yaml (preferred)
    parser.add_argument('-y', '--yaml', default=None, help='Name of configuration file. This is the preferred way to'
                                                           'pass parameters since it is easily reproducible. There are'
                                                           'also a lot of parameters')
    # load pickled data
    parser.add_argument('-l', '--load', action="store_true", help='Load stored pickle file (states.pl)')

    # args for everything in the yaml. These are ignored if yaml is passed
    # MD trajectory control
    parser.add_argument('-t', '--trajectory', default='PR_nojump.xtc', help='Path to input file.')
    parser.add_argument('-g', '--gro', default='em.gro', help='Name of .gro coordinate file.')
    parser.add_argument('-r', '--res', default='MET', help='Name of residue')
    parser.add_argument('-s', '--start', default=0, type=int, help='Frame at which to start analysis.')

    # association parameters
    parser.add_argument('-rc', '--coordinated_residue', default=None, help='Name of residue coordinated to residue_')
    parser.add_argument('-ac', '--coordinated_atoms', default=None, nargs='+', help='Name of residue coordinate to '
                        'residue')
    parser.add_argument('-ca', '--catoms', default=None, nargs='+', help='Name of atoms to calculate coordination '
                        'number with respect to. The center of mass will be used')
    parser.add_argument('-ta', '--atype', default=None, help='Element name of atoms of which you want coordination '
                                                             'number')
    parser.add_argument('-tc', '--coordinated_type', default=None, help='Element name of coordinated atoms')
    parser.add_argument('-cut', default=0.25, type=float, help='Maximum distance between pairs where they are considered'
                                                              'coordinated (nm)')

    # hbond parameters
    parser.add_argument('-hr', '--hresidues', nargs='+', default=['HII'], help='Residues to include in h-bond search')
    parser.add_argument('-ha', '--hatoms', action='append', nargs='+', help='Atoms to include for each '
                        'residue. Each list of atoms must be passed with a separate -a flag for each residue in'
                        'args.residues')
    parser.add_argument('-hd', '--hdistance_cut', default=.35, type=float, help='Maximum distance between acceptor and'
                                                                                'donor atoms')
    parser.add_argument('-hangle', '--hangle_cut', default=30, type=float, help='Maximum DHA angle to be considered an '
                                                                                'H-bond')
    parser.add_argument('-hacc', '--acceptors', default=False, nargs='+', help='If you only want hbonds with acceptor'
                        'atoms of a residue, use this flag. 1 for this condition, otherwise 0. Input as a list'
                        ' of 0 and 1 in the same order as args.residues')
    parser.add_argument('-hdonors', '--donors', default=False, nargs='+', help='If you only want hbonds with donor'
                        'atoms of a residue, use this flag. 1 for this condition, otherwise 0. Input as a list'
                        ' of 0 and 1 in the same order as args.residues')

    # partition parameters
    parser.add_argument('-ns', '--npts_spline', default=10, type=int, help='Number of points in each pore spline.')
    parser.add_argument('-rp', '--rpartition', default=0.75, type=float, help='Distance from pore center where system'
                                                                              'transitions from pore to tail region.')
    parser.add_argument('-ref', '--ref_atoms', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Names of atoms whose'
                        'centers of mass will be used to determine pore center.')

    return parser


class States:

    def __init__(self, traj, gro, res, start=0, **kwargs):
        """ Identify discrete states in a trajectory

        :param traj: name of GROMACS trajectory file (.xtc or .trr)
        :param gro: name of GROMACS coordinate file (.gro)
        :param res: name of residue whose states will be tracked
        :param start: frame to start analysis. Should be after system is equilibrated

        :type traj: str
        :type gro: str
        :type res: str
        :type start: int
        """

        self.residue = res

        # Read in all the parameters
        if 'association_params' not in kwargs:
            # default paramters for coordination_number.System. Can be modified with kwargs
            self.association_params = {'coordinated_residue': None, 'atoms': None, 'coordinated_atoms': None,
                                       'type': None, 'coordinated_type': None, 'begin': 0, 'end': -1, 'skip': 1,
                                       'com': True}
        else:
            self.association_params = kwargs['association_params']

        if 'hbond_params' not in kwargs:
            self.hbond_params = {'acceptors': False, 'residues': ['HII', 'MET'], 'donors': False, 'atoms': False,
                                 'angle_cut': 30, 'distance_cut': 0.35}
        else:
            self.hbond_params = kwargs['hbond_params']

        if 'partition_params' not in kwargs:
            self.partition_params = {'npts_spline': 10, 'r': 0.75, 'ref_atoms': ['C', 'C1', 'C2', 'C3', 'C4', 'C5']}
        else:
            self.partition_params = kwargs['partition_params']

        if 'state_labels' not in kwargs:
            self.state_labels = {'intails': 0, 'hbonded': 1, 'associated': 2, 'inpores': 3, 'state_dict': {'1000': 0,
                                 '1100': 1, '1010': 2, '1110': 3, '0001': 4, '0101': 5, '0011': 6, '0111': 7}}
        else:
            self.state_labels = kwargs['state_labels']

        if 'emission_params' not in kwargs:
            self.emission_params = {'distribution': 'Levy', 'lump_transitions': True}
        else:
            self.emission_params = kwargs['emission_params']

        self.nstates = len(self.state_labels['state_dict'].keys())

        print("Loading trajectory...", end='', flush=True)
        self.t = md.load(traj, top=gro)[start:]
        print("Done!")

        print("Identifying hydrogen bonds...")
        self.hbonds = self._identify_hydrogen_bonds(gro)
        print("Done!")

        print("Identifying other electrostatic associations...")
        self.associations = self._identify_associations(gro)
        print("Done!")

        print("Determining frame-by-frame partition of solutes between between pore and tail...")
        self.partition, self.com = self._define_partition()
        print("Done!")

        self.nT, self.nsolute = self.com.shape[:2]

        self.state_sequence = None
        self.transition_matrix = None
        self.count_matrix = None
        self.emissions = None
        self.fit_params = None
        self.hurst = None

        self._determine_state_sequence()
        self._make_transition_matrix()
        self.measure_state_emissions()
        self.calculate_hurst()

        self.box_length = self.t.unitcell_vectors[:, 2, 2].mean()

    def _identify_hydrogen_bonds(self, gro):
        """ Use hbonds.System to identify hydrogen bonds between the residue and the specified set of atoms

        :param gro: GROMACS coordinate file describing system

        :type gro: str
        """

        hb = hbonds.System(None, gro, t=self.t)

        residues = self.hbond_params['residues']
        if not self.hbond_params['acceptors']:
            acceptors = [False for _ in residues]
        else:
            acceptors = [bool(int(i)) for i in self.hbond_params['acceptors']]

        if not self.hbond_params['donors']:
            donors = [False for _ in residues]
        else:
            donors = [bool(int(i)) for i in self.hbond_params['donors']]

        if not self.hbond_params['atoms']:
            atoms = [['all'] for _ in residues]  # a default value
        else:
            atoms = self.hbond_params['atoms']

        while len(atoms) != len(residues):
            atoms.append(['all'])

        for i, r in enumerate(residues):
            hb.set_eligible(r, atoms[i], acceptor_only=acceptors[i], donors_only=donors[i])

        hb.identify_hbonds(self.hbond_params['distance_cut'], self.hbond_params['angle_cut'])

        # simplify list to something binary. True indicates that the residue is involved in an hbond in some capacity
        res_number, nres = hb.number_residues(self.residue)
        nT = len(hb.hbonds)
        hb_binary = np.zeros([nT, nres], dtype=bool)  # boolean to save memory
        keys = set(res_number.keys())
        for t in range(self.t.n_frames):
            bonds = set(hb.hbonds[t][:3, :].flatten())  # find which hbonding atoms are part of res
            vals = np.unique([res_number[int(x)] for x in keys.intersection(bonds)])  # identify unique res hbonds
            try:
                hb_binary[t, vals] = True
            except IndexError:  # some frames might not have any hbonding
                pass

        return hb_binary

    def _identify_associations(self, gro):
        """ Use coordination_number.System to identify electrostatic associations between the residue and the specified
        set of atoms

        :param gro: GROMACS coordinate file describing system

        :type gro: str
        """

        association = coordination_number.System(None, gro, residue=self.residue,
                                                 coordinated_residue=self.association_params['coordinated_residue'],
                                                 atoms=self.association_params['atoms'],
                                                 coordinated_atoms=self.association_params['coordinated_atoms'],
                                                 type=self.association_params['type'],
                                                 ctype=self.association_params['coordinated_type'],
                                                 com=self.association_params['com'], t=self.t)

        # calculate pairwise distance between all points in self.com and self.com_coordinated
        association.distance_search(cut=self.association_params['cut'])

        association.n_coordinated(plot=False)

        return association.ncoord.astype(bool)

    def _define_partition(self):
        """ Classify solutes as either {True: in the pores} or {False: in the tails} based on their center of mass
        positions

        :param npts: number of points in spline

        :type npts: int

        :return: array of True/False in/out of pore region
        """

        com = self._center_of_mass()
        atoms = [a.index for a in self.t.topology.atoms if a.name in self.partition_params['ref_atoms']]
        pore_centers = physical.trace_pores(self.t.xyz[:, atoms, :], self.t.unitcell_vectors,
                                            self.partition_params['npts_spline'], npores=4, progress=True, save=True,
                                            savename='spline.pl')[0]
        return physical.partition(com, pore_centers, self.partition_params['r'], unitcell=self.t.unitcell_vectors,
                                  npores=4, spline=True), com

    def _center_of_mass(self):

        print('Calculating center of mass trajectories of residue %s' % self.residue)

        residue = topology.Residue(self.residue)  # get resiude attributes

        ndx = [a.index for a in self.t.topology.atoms if a.residue.name == self.residue]  # index of all residue atoms
        names = [a.name for a in self.t.topology.atoms if a.residue.name == res][
                :residue.natoms]  # names of atoms in one residue
        mass = [residue.mass[x] for x in names]  # mass of atoms in order that they appear in file

        return physical.center_of_mass(self.t.xyz[:, ndx, :], mass)  # determine center of mass trajectories

    def _translate_state(self, s):
        """ convert list of True/False to string of 1's and 0's and use that to determine its state label from
        self.state_labels['state_dict']

        :param s: tuple of True/False
        :type s: tuple

        :return: state
        :rtype: int
        """

        return self.state_labels['state_dict'][''.join([str(int(i)) for i in s])]

    def _determine_state_sequence(self):
        """ Based on trapping mechanisms, assign a state to each time point and each solute
        """

        print("Determining state sequence...")
        self.state_sequence = np.zeros_like(self.hbonds, dtype=int)

        # note that these are in same order as zip below. This implies that self.partition_params['labels'] also has
        # a specific order in the .yaml
        labels = [self.hbond_params['label'], self.association_params['label'], self.partition_params['labels'][0],
                  self.partition_params['labels'][1]]

        order = np.argsort([self.state_labels[i] for i in labels])

        for t in tqdm.tqdm(range(self.nT)):

            for i, s in enumerate(zip(self.hbonds[t], self.associations[t], ~self.partition[t], self.partition[t])):
                self.state_sequence[t, i] = self._translate_state(np.array(s)[order])

        print("Done!")

    def _make_transition_matrix(self, start=1, step=1):
        """ Estimate a transition matrix based on the state sequence
        """

        self.transition_matrix = np.zeros([self.nstates, self.nstates])
        self.count_matrix = np.zeros_like(self.transition_matrix, dtype=int)

        for t in range(start, self.nT, step):  # start at frame 1. May need to truncate more as equilibration
            transitioned_from = self.state_sequence[t - 1, :]
            transitioned_to = self.state_sequence[t, :]
            for pair in zip(transitioned_from, transitioned_to):
                self.count_matrix[pair[0], pair[1]] += 1

        # normalize so rows sum to unity
        self.transition_matrix = (self.count_matrix.T / self.count_matrix.sum(axis=1)).T

    def verify_markovinity(self, timesteps=(1, 2, 4), cmap='bwr', dt=0.5, show=False, save=False):
        """ Verify that the transition matrix is Markovian

        The condition of detailed balance is satisfied if the following holds:
        .. math::

            T_{i,j}P_j(t=\infty) = T_{j, i}P_j(t=\infty)

        In other words, the off-diagonal entries of the transition matrix are equal.

        A further test is to show that the transition matrix is independent of the time step.

        :param timesteps: timesteps at which to create transition matrix. A timestep of x would add to the transition
        matrix every x frames.

        :type timesteps: tuple
        """

        start = 0.5
        end = self.nstates + start
        extent = [start, end, start, end]
        x_positions = np.linspace(start, end - 1, self.nstates)
        y_positions = np.linspace(start, end - 1, self.nstates)

        fig = plt.figure(figsize=(15, 6))
        fig2 = plt.figure(figsize=(15, 6))
        grid = AxesGrid(fig, 111, nrows_ncols=(1, 3), cbar_mode='none', cbar_location='top', cbar_pad=0.2)
        grid2 = AxesGrid(fig2, 111, nrows_ncols=(1, 3), cbar_mode='none', cbar_location='top', cbar_pad=0.2)

        for s, ax in enumerate(grid):

            self._make_transition_matrix(step=timesteps[s])

            T = self.transition_matrix[::-1, :].T
            C = self.count_matrix[::-1, :].T

            ax.imshow(T.T, extent=extent, origin='lower', cmap=cmap)
            grid2.axes_all[s].imshow(C.T, extent=extent, origin='lower', cmap=cmap,
                                     norm=colors.LogNorm(self.count_matrix.min(), self.count_matrix.max()))

            for yndx, y in enumerate(y_positions):
                for xndx, x in enumerate(x_positions):

                    label = T[xndx, yndx]
                    counts = C[xndx, yndx]
                    ax.text(x + start, y + start, '%.2f' % label, color='black', ha='center', va='center', weight='bold')
                    grid2.axes_all[s].text(x + start, y + start, '%d' % counts, color='black', ha='center', va='center', weight='bold')

            ax.tick_params(labelsize=14)
            ax.set_xticks(np.arange(1, self.nstates + 1))
            ax.set_title('Timestep = %.1f ns' % (dt * timesteps[s]), fontsize=14)
            ax.set_yticklabels(np.arange(9, 0, -1))

            grid2.axes_all[s].tick_params(labelsize=14)
            grid2.axes_all[s].set_xticks(np.arange(1, self.nstates + 1))
            grid2.axes_all[s].set_title('Timestep = %.1f ns' % (dt * timesteps[s]), fontsize=14)
            grid2.axes_all[s].set_yticklabels(np.arange(9, 0, -1))

        fig.tight_layout()
        fig2.tight_layout()

        if save:
            fig.savefig('%s_transitions.pdf' % self.residue)
            fig2.savefig('%s_counts.pdf' % self.residue)
        if show:
            plt.show()

    def measure_state_emissions(self, dim=2, fit_function=levy):
        """ Measure observations as a function of state.

        :param dim: dimension of trajectory whose emissions will be measured
        :param fit_function: name of function to fit emissions to. Only works for Levy currently due to naming

        :type dim: int or list
        :type fit_function: scipy.stats.rv_continuous
        """

        self.emissions = [[] for _ in range(self.nstates)]  # row for each state
        for row in range(self.nstates):  # column for each state
            self.emissions[row] = [[] for _ in range(self.nstates)]

        print("Measuring state emissions...", flush=True, end='')
        for t in range(1, self.nT):
            for e in range(self.nsolute):
                self.emissions[self.state_sequence[t - 1, e]][self.state_sequence[t, e]].append(self.com[t, e, dim] -
                                                                                                self.com[t - 1, e, dim])
        print('Done!')

        if self.emission_params['lump_transitions']:
            self._lump_transition_emissions()

        self._fit_emissions_distributions(fit_function)

    def _lump_transition_emissions(self):
        """ Lump all of the emission distributions for transitions between states together. The last entry in the
        emissions list
        """

        emissions = [[] for _ in range(self.nstates + 1)]  # a distribution for each state and 1 for transitions
        for s1 in range(self.nstates):
            for s2 in range(self.nstates):
                if s1 == s2:
                    emissions[s1] += self.emissions[s1][s2]
                else:
                    emissions[-1] += self.emissions[s1][s2]

        self.emissions = emissions

    def _fit_emissions_distributions(self, fit_function, nboot=200):
        """ Fit Levy distribution to emission distributions

        :param fit_function: functional form of emission distributions
        :param nboot: number of bootstrap trials
        """

        print('Fitting emission distributions...', flush=True, end='')
        if self.emission_params['lump_transitions']:
            self.fit_params = np.zeros([nboot, self.nstates + 1, 3])
            for s in tqdm.tqdm(range(self.nstates + 1), unit='States'):
                for b in range(nboot):
                    emissions = np.random.choice(self.emissions[s], size=len(self.emissions[s]), replace=True)
                    self.fit_params[b, s, :] = fit_function.fit_levy(emissions, beta=0)[0].x  # alpha, mu, sigma
        else:
            self.fit_params = np.zeros([nboot, self.nstates, self.nstates, 3])  # only 3 parameters since beta is fixed

            for i in tqdm.tqdm(range(self.nstates)):
                for j in range(self.nstates):
                    for b in range(nboot):
                        emissions = np.random.choice(self.emissions[i][j], size=len(self.emissions[i][j]), replace=True)
                        self.fit_params[b, i, j, :] = fit_function.fit_levy(emissions, beta=0)[0].x

        print('Done!')

    def maximum_emission(self):
        """ Find the largest observed emission

        :return The largest observed emission
        :rtype: float
        """

        max_emissions = []
        if self.emission_params['lump_transitions']:
            for i in self.emissions:
                max_emissions.append(max(np.abs(i)))
        else:
            for i in self.emissions:
                for j in i:
                    max_emissions.append(max(np.abs(j)))

        return max(max_emissions)

    def plot_emissions(self, save=False, savename='emissions'):
        """ Plot emissions within discrete states. (This does not plot transition distributions)

        :param save: save plot when done
        :param savename: name to give plot

        :type save: bool
        :type savename: str
        """

        fig, ax = plt.subplots(2, 4, figsize=(12, 7))  # specific to 8 states
        mechanisms = list(self.state_labels['state_dict'].keys())
        states = [i for i in self.state_labels.keys() if i != 'state_dict']
        map_state_labels = {self.state_labels[i]: i for i in states}
        x = np.linspace(-1, 1, 1000)

        for s in range(self.nstates):

            ax[s // 4, s % 4].hist(self.emissions[s][s], bins=50, density=True)
            string_mechanism = '+'.join([map_state_labels[i] for i, x in enumerate(list(mechanisms[s])) if x != '0'])
            ax[s // 4, s % 4].set_title('State %d\n %s\n no. data points: %d' % ((s + 1), string_mechanism,
                                                                                 self.count_matrix[s, s]))  # count matrix since only self-transitions are used
            ax[s // 4, s % 4].set_xlim(-1, 1)
            ax[s // 4, s % 4].plot(x, levy_stable.pdf(x, self.fit_params[s, s, 0], 0,
                                   loc=self.fit_params[s, s, 1], scale=self.fit_params[s, s, 2]), linewidth=2)
            ax[s // 4, s % 4].tick_params(labelsize=14)

        plt.tight_layout()

        if save:
            plt.savefig('%s.pdf' % savename)

        plt.show()

    def transition_autocorrelation(self, nboot=200, confidence=95, plot=False, max_k=15, fontsize=14):
        """ Calculate the autocorrelation function of hops only when transitions occur.

        :param nboot: number of bootstrap trials for generating confidence intervals
        :param confidence: percent confidence interval (out of 100)
        :param plot: plot autocorrelation function
        :param fontsize: size of font on plot axes
        :param max_k: largest time lag to plot

        :type nboot: int
        :type confidence: float
        :type fontsize: int
        :type plot: bool
        :type max_k: int

        :return: acf
        :rtype: numpy.ndarray
        """

        transitions = []
        for s in range(self.nsolute):
            switch_points = timeseries.switch_points(self.state_sequence[:, 0])
            transitions.append(self.com[switch_points[:-1] + 1, s, 2] - self.com[switch_points[:-1], s, 2])

        if plot:

            # calculate the autocorrelation function of uneven length trajectories
            acf, errorbars = timeseries.acf_uneven(transitions, nboot=nboot, confidence=confidence)
            timeseries.plot_autocorrelation(acf, errorbars=errorbars, bootstrap=False, show=True, max_k=max_k,
                                            fontsize=fontsize)

        else:

            acf = timeseries.acf_uneven(transitions)  # get acf of each traj without bootstrapping

        return acf

    def state_autocorrelation(self):
        """ Calculate autocorrelation while in a specific state
        """

        sequences = [[] for _ in range(self.nstates)]
        for solute in range(self.nsolute):
            ndx = 0
            for state, group in itertools.groupby(self.state_sequence[:, solute]):
                n = len(list(group))  # sequence length
                positions = self.com[ndx:(ndx + n), solute, 2]  # z-positions during dwell sequence
                if positions.size > 1:  # need at least two data points to get a hop within a state
                    sequences[state].append(positions[1:] - positions[:-1])  # sequence of increments during dwell
                ndx += n

        acfs = []

        for seq in sequences:
            acfs.append(timeseries.acf_uneven(seq))

        # for i, seq in enumerate(sequences):
        #     acf, errorbars = timeseries.acf_uneven(seq, nboot=200)
        #     timeseries.plot_autocorrelation(acf, show=False, errorbars=errorbars)
        #     plt.figure()
        #     biggest = np.argmax([len(j) for j in sequences[i]])
        #     plt.plot(sequences[i][biggest])
        #     plt.figure()
        #     plt.plot(np.cumsum(sequences[i][biggest]))
        #     plt.show()
        # exit()

        return acfs

    def calculate_hurst(self, nboot=200, max_k=15, plot=False, state=0):
        """
        """

        self.hurst = np.zeros([self.fit_params.shape[0], nboot])

        for h, acf in enumerate(self.state_autocorrelation()):
            self.hurst[h, :] = timeseries.hurst(acf, nboot=nboot, max_k=max_k)

        acf = self.transition_autocorrelation()
        self.hurst[-1] = timeseries.hurst(acf, nboot=nboot, max_k=max_k)

        if plot:

            acf = self.state_autocorrelation()[state]
            weights = np.array([np.nonzero(acf[i, :])[0].size for i in range(acf.shape[0])])
            fig, ax = plt.subplots(1, 1)
            ax.plot(np.average(acf, axis=0, weights=weights)[:max_k], linewidth=2, label='Simulated autocorrelation')
            ax.plot(np.arange(max_k), fitting_functions.hurst_autocovariance(np.arange(max_k),
                    np.mean(self.hurst[state])), '--', color='black', lw=2, label='Fit theoretical autocorrelation')
            ax.legend(fontsize=14)

            plt.show()


class Chain:

    def __init__(self, count_matrix, emission_parameters, hurst_parameters=None, emission_function=levy_stable,
                 bound=7.6, speedup=True, fbm=False):
        """ Generate markov chains based on a probability transition matrix and some emission paramters describing the
        observations.

        :param count_matrix: n x n matrix of counts of transitions between each of the n states. This will be
        converted to a transition matrix.
        :param emission_parameters: parameters describing distribution of emissions for each state (n + 1 x p) where p
        is the number of paramters, and n the number of states. The extra set of parameters should be those of the
        lumped together transition emission distribution.
        :param hurst_parameters: distribution of self-similarity parameters for each state and transitions between
        states. (nstates + 1, n) where n is the number of hurst parameters in the distribution. Each distribution is
        created by bootstrapping fits to the acf.
        :param emission_function: name of probability distribution describing emissions

        :type transition_matrix: np.ndarray
        :type emission_parameters: np.ndarray
        :type hurst_parameters: np.ndarray
        :type emission_function: scipy.stats.continuous_rv
        """

        self.count_matrix = count_matrix
        self.transition_matrix = (self.count_matrix.T / self.count_matrix.sum(axis=1)).T
        self.emission_parameters = emission_parameters
        self.emission_function = emission_function
        self.hurst_parameters = hurst_parameters

        self.nstates = self.transition_matrix.shape[0]

        # self.lumped_transitions = False
        # if self.emission_parameters.shape[0] == (self.transition_matrix.shape[0] + 1):
        #     self.lumped_transitions = True

        self.chains = None
        self.states = None
        self.msds = None
        self.limits = None
        self.bound = bound
        self.timeout = None
        self.fbm = fbm
        self.gaussian_params = None

        self.average_params = None
        self.truncated_distributions = None
        self.flms = None

        if self.hurst_parameters is None:

            raise Exception('You must provide a distribution of Hurst parameters')

        elif self.fbm:

            self._convert_levy_to_gaussian()

        else:

            print('Generating Interpolator Objects ...', end='', flush=True)
            self.hurst_interpolator = flm_sim_params.HurstCorrection()
            self.truncation_interpolator = flm_sim_params.TruncateLevy()
            print('Done!')

            self.m = 256
            self.Mlowerbound = 6000

            if speedup:

                self._get_average_params()

                self._get_truncated_distributions()

    def _convert_levy_to_gaussian(self):
        """ Convert truncated levy distribution to a Gaussian distribution by taking the standard deviation of the
        truncated distribution. Assumes zero skew and zero mean.
        """

        self.gaussian_params = np.zeros(self.nstates + 1)

        print('Converting truncated Levy distribution to Gaussian...', flush=True)
        for s in tqdm.tqdm(range(self.nstates + 1)):

            alpha, scale = self.emission_parameters[s, [0, 2]]

            levy = sampling.TruncatedLevyDistribution(self.bound, alpha, scale, 1000000, beta=0, mean=0)

            self.gaussian_params[s] = np.std(levy.empirical_distribution)

        print('Done!')

    def _get_average_params(self):

        print('Interpolating Mean Parameters ...', flush=True, end='')

        self.average_params = {'alpha': {}, 'scale': {}, 'hurst': {}, 'corrected_bound': {}, 'corrected_H': {}}

        # This loop could take some time
        for s in tqdm.tqdm(range(self.nstates + 1)):

            H = np.mean(self.hurst_parameters[s])
            alpha, scale = self.emission_parameters[s, [0, 2]]
            corrected_bound = self.truncation_interpolator.interpolate(H, alpha, self.bound, scale)
            corrected_H = self.hurst_interpolator.interpolate(H, alpha)

            self.average_params['hurst'][s] = H
            self.average_params['alpha'][s] = alpha
            self.average_params['scale'][s] = scale
            self.average_params['corrected_bound'][s] = corrected_bound
            self.average_params['corrected_H'][s] = corrected_H

        print('Done!')

    def _get_truncated_distributions(self):

        # BJC: I think this might be wrong. You should get alpha and scale from the original emission params

        print('Creating empirical truncated Levy stable distributions ...', flush=True, end='')
        self.truncated_distributions = []
        for s in tqdm.tqdm(range(self.nstates + 1)):
            self.truncated_distributions.append(
                sampling.TruncatedLevyDistribution(self.bound, self.average_params['alpha'][s],
                                                   self.average_params['scale'][s], 1000000, beta=0, mean=0))
            # self.average_params['scale'][s] = np.std(self.truncated_distributions[s].sample(10000))  # REMOVE!
        print('Done!')

    def _initialize_flm(self):

        print('Initializing FLM objects for each state ...', flush=True, end='')
        self.flms = []
        for s in tqdm.tqdm(range(self.nstates + 1)):

            corrected_H = self.average_params['corrected_H'][s]
            alpha = self.average_params['alpha'][s]
            scale = self.average_params['scale'][s]
            corrected_bound = self.average_params['corrected_bound'][s]
            #
            # self.flms.append(FLM(corrected_H, alpha, m=self.m, M=4, N=l, scale=scale, truncate=corrected_bound,
            # correct_hurst=False, correct_truncation=False)

    def generate_realizations(self, n, length, bound=7.6, m=256, Mlowerbound=6000, nt=1, speedup=True, quiet=False,
                              timeout=None):
        """ Generate Markov chains using transition matrix and emission probabilities

        :param n: number of independent trajectories to generate
        :param length: length (number of steps) of generated trajectories
        :param bound: largest hop allowable. This essentially truncates the emission distribution
        :param m: mesh size. Choose a value that is a power of 2
        :param Mlowerbound: M is the kernel function cut-off. This is the lowest possible value of M that will be \
        chosen. M will be chosen such that M + length is a power of 2. Higher values of M will lead to more accurate \
        calculation of the correlation structure at large time lags.
        :param nt: number of threads to use for trajectory realization generation.
        :param speedup: speed up trajectory generation by using the mean of the parameter distributions so that
        interpolation doesn't need to be done repeatedly.

        :type n: int
        :type length: int
        :type bound: float
        :type m: int
        :type Mlowerbound: int
        :type nt: int
        :type speedup: bool
        """

        # allow resetting of bound from that passed to __init__
        if np.abs(self.bound - bound) > 0.0001 and not self.fbm:

            self.bound = bound

            if speedup:

                self._get_average_params()

                self._get_truncated_distributions()

        # if bound isn't reset but speedup is passed and it wasn't passed to __init__
        elif speedup and self.average_params is None and not self.fbm:

            self._get_average_params()

            self._get_truncated_distributions()

        self.timeout = timeout

        trajectories_per_thread = n // nt

        self.chains = np.zeros([length, trajectories_per_thread * nt])
        self.states = np.zeros([length, trajectories_per_thread * nt], dtype=int)

        self.m = m
        self.Mlowerbound = Mlowerbound

        if not quiet:
            print("Generating trajectories...")

        if nt > 1:

            pool = Pool(processes=nt)  # set up multiprocessing pool

            trajectories_per_thread = n // nt

            arguments = [(trajectories_per_thread, length, bound, quiet) for _ in range(nt)]  # arguments for each thread

            result = pool.starmap(self._realizations, arguments)  # do the calculations

            # unpack result list
            for i, r in enumerate(result):

                self.chains[:, i * trajectories_per_thread:(i + 1) * trajectories_per_thread] = r[0]
                self.states[:, i * trajectories_per_thread:(i + 1) * trajectories_per_thread] = r[1]

            pool.close()
            pool.join()

        else:

            self.chains, self.states = self._realizations(n, length, bound, quiet)

    def _realizations(self, n, length, bound, quiet):
        """ Generate multiple realization of Markov chain. This function exists by itself to aid parallelization

        :param n: number of trajectories to generate
        :param length: length of trajectory
        :param bound: largest hop allowable. This essentially truncates the emission distribution

        :type n: int
        :type length: int
        :type bound: float

        :return: chains, states
        """

        np.random.seed()

        chains = np.zeros([length, n])
        states = np.zeros([length, n], dtype=int)
        for i in tqdm.tqdm(range(n), disable=quiet):
            chains[:, i], states[:, i] = self._realization(length, bound=bound)

        return chains, states

    def _realization(self, length, bound=7.6):
        """ Generate single realization of Markov chain

        :param length: length of trajectory
        :param bound: largest hop allowable. This essentially truncates the emission distribution

        :type length: int
        :type bound: float
        """

        if self.hurst_parameters is None:

            return self._uncorrelated_realization(length, bound=bound)

        else:

            return self._correlated_realization(length, bound=bound)

    def _correlated_realization(self, length, bound=7.6):
        """ Generate a realization with correlated emissions

        :param length: length of trajectory
        :param bound: largest hop allowable. This essentially truncates the emission distribution

        :type length: int
        :type bound: float
        """

        traj = np.zeros([length])
        states = np.zeros([length], dtype=int)

        # Determine state sequence
        states[0] = np.random.randint(self.nstates)  # choose initial state with uniform probability

        for t in range(1, length):

            p = np.random.dirichlet(self.count_matrix[states[t], :])
            #states[t] = np.random.choice(self.nstates, p=self.transition_matrix[states[t - 1], :])
            states[t] = np.random.choice(self.nstates, p=p)

        # find where transitions occur
        transitions = timeseries.switch_points(states)  # indices of state array where transitions occur
        if states[0] == states[1] and transitions[0] == 0:  # switch_points always includes first and last data point
            transitions = transitions[1:]

        for i, t in enumerate(transitions[:-1]):

            dwell = transitions[i + 1] - t  # plus two since using enumerate and starting at index 1 of transitions
            if dwell > 1:
                if dwell > 2:
                    traj[t + 1: transitions[i + 1]] = self._flm_sequence(dwell - 1, states[t], bound=bound)
                elif dwell == 2:

                    # one data point so generate single point from emission distribution to save time
                    if self.fbm:  # this section could definitely be written better

                        rv = np.random.normal(0, self.gaussian_params[states[t]], 1)

                    else:

                        rv = bound + 1
                        while np.abs(rv) > bound:

                                rv = self.emission_function.rvs(self.emission_parameters[states[t], 0], 0,
                                                                loc=self.emission_parameters[states[t], 1],
                                                                scale=self.emission_parameters[states[t], 2])

                    traj[t + 1] = rv
                    #traj[t + 1] = self._flm_sequence(1, states[t], bound=bound)

        # slow step
        traj[transitions[:-1]] = self._flm_sequence(transitions.size - 1, self.nstates, bound=bound)

        return np.cumsum(traj), states

    def _flm_sequence(self, l, state, bound=7.6):
        """ Generate a sequence of correlated draws from a levy process

        :param l: length of sequence
        :param state: state
        :param bound: largest allowable hop

        :type l: int
        :type state: int
        :type bound: float

        :return: sequence
        :rtype: numpy.ndarray
        """

        if not self.fbm:

            if state != 8:

                return np.zeros(l)

            else:

                # for efficiency, ensure m * (l + M) is a power of 2
                M = l + self.Mlowerbound
                n = np.log(M) / np.log(2)
                M = int(2 ** np.ceil(n) - l)

                if self.average_params is not None:

                    alpha = self.average_params['alpha'][state]
                    scale = self.average_params['scale'][state]
                    corrected_H = self.average_params['corrected_H'][state]
                    corrected_bound = self.average_params['corrected_bound'][state]

                else:

                    H = np.random.choice(self.hurst_parameters[state])  # random hurst parameter from transition distribution
                    alpha, scale = self.emission_parameters[state, [0, 2]]
                    corrected_bound = self.truncation_interpolator.interpolate(H, alpha, bound, scale)
                    corrected_H = self.hurst_interpolator.interpolate(H, alpha)

                trunc_dist = None
                if self.truncated_distributions is not None:
                    trunc_dist = self.truncated_distributions[state]

                flm = FLM(corrected_H, alpha, m=self.m, M=M, N=l, scale=scale, truncate=corrected_bound,
                          correct_hurst=False, correct_truncation=False)

                # generate a single realization
                flm.generate_realizations(1, progress=False, truncated_distribution=trunc_dist)

                #H = self.average_params['hurst'][state]

                #return scale * fbm.FBM(int(l), H, method="daviesharte").fgn() / ((1.0 / l) ** H)

                return flm.noise[0, :l]

        else:

            H = np.random.choice(self.hurst_parameters[state])
            scale = self.gaussian_params[state]

            if state != 8:

                return np.zeros(l)

            else:

                return scale * fbm.FBM(int(l), H, method="daviesharte").fgn() / ((1.0 / l) ** H)

        # if state == 8:
        #     # flm = FLM(corrected_H, alpha, m=self.m, M=M, N=l, scale=scale, truncate=corrected_bound,
        #     #           correct_hurst=False, correct_truncation=False)
        #     #
        #     # # generate a single realization
        #     # flm.generate_realizations(1, progress=False, truncated_distribution=trunc_dist)
        #     # #(sigma, alpha, hurst, max_hop)
        #     H = 0.35
        #
        #     return scale * fbm.FBM(l, H, method="daviesharte").fgn() / ((1.0 / l) ** H)
        #
        #     #return flm.noise[0, :l]
        #
        # else:
        #
        #     return np.zeros(l)

        # flm.autocorrelation()
        # H = np.log(2 * flm.acf[:, 1].mean() + 2) / (2 * np.log(2))
        # print(H)
        # plt.hist(flm.noise[0, :], bins=50, range=(-bound/2, bound/2), density=True)
        # x = np.linspace(-bound/2, bound/2, 200)
        # plt.plot(x, levy_stable.pdf(x, alpha=alpha, loc=0, beta=0, scale=scale))
        # plt.show()
        # exit()

    def _uncorrelated_realization(self, length, bound=7.6):
        """ Generate a trajectory with uncorrelated emissions

        :param length: length of trajectory
        :param bound: largest hop allowable. This essentially truncates the emission distribution

        :type length: int
        :type bound: float
        """

        traj = np.zeros([length])
        states = np.zeros([length], dtype=int)

        states[0] = np.random.randint(self.nstates)  # choose initial state with uniform probability
        hop = self.hop(states[0], states[0])  # assume hop drawn from emission distribution in that state
        while abs(hop) > bound:
            hop = self.hop(states[0], states[0])
        traj[0] = hop

        for t in range(length - 1):

            hop = self.hop(states[t - 1], states[t])
            while abs(hop) > bound:
                hop = self.hop(states[t - 1], states[t])
            traj[t] = hop

            states[t + 1] = np.random.choice(self.nstates, p=self.transition_matrix[states[t], :])

        hop = self.hop(states[-2], states[-1])
        while abs(hop) > bound:
            hop = self.hop(states[-2], states[-1])
        traj[-1] = hop

        return np.cumsum(traj), states

    def hop(self, state1, state2):
        """ Draw a hop from the appropriate emission distribution depending on the current and previous state

        :param state1: previous state
        :param state2: new state

        :type state1: int
        :type state2: int
        """

        if state1 == state2:

            return self.emission_function.rvs(self.emission_parameters[state1, 0], 0,
                                              loc=self.emission_parameters[state1, 1],
                                              scale=self.emission_parameters[state1, 2])
        else:

            return self.emission_function.rvs(self.emission_parameters[-1, 0], 0,
                                              loc=self.emission_parameters[-1, 1],
                                              scale=self.emission_parameters[-1, 2])

    def calculate_msd(self, nboot=200):
        """ Calculate the average mean squared displacement (MSD) of simulated particle trajectories

        :param nboot: number of bootstrap trials for generating statistics

        :type nboot: int
        """

        self.msds = timeseries.msd(self.chains[..., np.newaxis], 0)
        self.limits = timeseries.bootstrap_msd(self.msds, nboot)

    def verify_transition_matrix(self, start=1):
        """ Re-calculate the transition matrix based on the generated state sequences. This is meant to check that you
        are properly sampling the real transition matrix.
        """

        transition_matrix = np.zeros([self.nstates, self.nstates])
        count_matrix = np.zeros_like(transition_matrix, dtype=int)

        for t in range(start, self.states.shape[0]):  # start at frame 1. May need to truncate more as equilibration
            transitioned_from = self.states[t - 1, :]
            transitioned_to = self.states[t, :]
            for pair in zip(transitioned_from, transitioned_to):
                count_matrix[pair[0], pair[1]] += 1

        # normalize so rows sum to unity
        transition_matrix = (count_matrix.T / count_matrix.sum(axis=1)).T

        print(transition_matrix)
        print(self.transition_matrix)

    def plot_msd(self, cutoff=0.4, dt=0.5, save=False, savename='msd.pdf', label='MSD', overlay=False, show=True):
        """ Plot mean squared displacement curve.

        :param cutoff: cut off the plot at this fraction of the total trajectory
        :param dt: time step (ns)
        :param save: save plot
        :param savename: name of saved plot
        :param label: legend label for plot
        :param overlay: if False, create a new figure
        :param show: show plot after

        :type cutoff: float
        :type dt: float
        :type save: bool
        :type savename: str
        :type label: str
        :type overlay: bool
        :type show: bool
        """

        if not overlay:
            plt.figure()

        time = np.arange(self.msds.shape[0])*dt
        end = int(len(time)*cutoff)
        avg = self.msds.mean(axis=1)

        plt.plot(time[:end], avg[:end], label=label, linewidth=2)
        plt.fill_between(time[:end], avg[:end] + self.limits[0, :end], avg[:end] - self.limits[1, :end], alpha=0.7)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('MSD (nm$^2$)', fontsize=14)

        print("MSD after %.1f ns: %.1f CI (%.2f, %.2f)" % (time[end], avg[end], avg[end] - self.limits[1, end],
                                                           avg[end] + self.limits[0, end]))
        plt.tight_layout()
        if save:
            plt.savefig(savename)
        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if not args.load:

        if not args.yaml:

            association_params = {'coordinated_residue': args.coordinated_residue, 'atoms': args.catoms,
                                  'coordinated_atoms': args.coordinated_atoms, 'type': args.type,
                                  'coordinated_type': args.coordinated_type,
                                  'begin': 0, 'end': -1, 'skip': 1, 'com': True}

            hbond_params = {'acceptors': args.acceptors, 'residues': args.hresidues, 'donors': args.donors,
                            'atoms': args.hatoms, 'angle_cut': args.hangle_cut, 'distance_cut': args.hdistance_cut}

            partition_params = {'npts_spline': args.npts_spline, 'r': args.rpartition, 'ref_atoms': args.ref_atoms}

            state_labels = {'intails': 0, 'hbonded': 1, 'associated': 2, 'inpores': 3, 'state_dict': {'1000': 0,
                            '1100': 1, '1010': 2, '1110': 3, '0001': 4, '0101': 5, '0011': 6, '0111': 7}}

            traj, gro, res, start = args.trajectory, args.gro, args.res, args.start

        else:  # preferred

            with open(args.yaml, 'r') as yml:
                cfg = yaml.load(yml)

            trajectory = cfg['trajectory']
            traj, gro, res, start = trajectory['trajectory'], trajectory['gro'], trajectory['res'], trajectory['start']

            association_params = cfg['association_params']
            hbond_params = cfg['hbond_params']
            partition_params = cfg['partition_params']
            state_labels = cfg['state_labels']
            emission_params = cfg['emission_params']

            if emission_params['distribution'].lower() != 'levy':
                sys.exit("The emission distribution, '%s', is not implemented" % emission_params['distribution'].lower())

        states = States(traj, gro, res, start=start, association_params=association_params, hbond_params=hbond_params,
                        partition_params=partition_params, state_labels=state_labels, emission_params=emission_params)

        states.t = None  # to save memory
        file_rw.save_object(states, 'states.pl')

        #states.plot_emissions()

    else:

        print("Loading pickled object...", flush=True, end='')
        states = file_rw.load_object('states.pl')
        print("Done!")

        #states.measure_state_emissions(fit_function=levy)
        #states.plot_emissions()

    # fig, ax = plt.subplots(8, 8, figsize=(10, 10), sharex=True, sharey=True)
    # for i in range(states.nstates):
    #     for j in range(states.nstates):
    #         ax[i, j].hist(states.emissions[i][j], bins=50, density=True, label='n = %d' % len(states.emissions[i][j]))
    #         ax[i, j].legend(fontsize=8)
    #         ax[i, j].tick_params(labelsize=8)
    #
    # plt.xlim(-0.5, 0.5)
    # plt.tight_layout()
    # plt.show()
    # exit()
    # print(states.state_sequence)
    # exit()

    states.verify_markovinity()

    #states._make_transition_matrix(1, end=100)

    total = states.count_matrix.sum(axis=0)
    p = total / total.sum()  # probability of being in state i at any time t
    w, v = np.linalg.eig(states.transition_matrix.T)

    print(w[0])
    print(v[:, 0] / v[:, 0].sum())
    print(p)
    plt.show()
    exit()
    ############ Cut-off Justification Plot ##################
    # alpha, mu, sigma = states.fit_params[-1]
    # bins, edges = np.histogram(states.emissions[-1], bins=200, range=(-1.75, 1.75), density=True)
    # bin_width = edges[1] - edges[0]
    # centers = edges[:-1] + bin_width / 2
    # ratio = [bins[i] / levy_stable.pdf(x, alpha=alpha, beta=0, loc=mu, scale=sigma) for i, x in enumerate(centers)]
    #
    # plt.figure(figsize=(10, 6))
    # plt.bar(centers, bins, bin_width, align='center', alpha=0.5, label='Empirical Disribution', color='xkcd:blue')
    # plt.plot(centers, levy_stable.pdf(centers, alpha=alpha, beta=0, loc=mu, scale=sigma), '--', lw=2,
    #          color='xkcd:orange', label='Analytical PDF')
    # plt.plot(centers, ratio, label='Ratio of Empirical:Analytical PDF', color='red', lw=2)
    # plt.plot(np.linspace(-1.75, 1.75, 3), np.ones([3]), '--', color='black', label='Exact agreement')
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    # plt.legend(fontsize=14)
    # plt.ylabel('Probability Density', fontsize=14)
    # plt.tight_layout()
    # plt.show()
    # exit()
    ############################################################
    ###### Show that Levy is better than Gaussian ######
    ###### Numbers that show undersampling of tails #######

    # alpha, mu, sigma = states.fit_params[-1]
    # transition_emissions = np.array(states.emissions[-1])
    # print(transition_emissions.size)
    # print(np.where(np.abs(transition_emissions) > 0.8)[0].size)
    #
    # plt.hist(states.emissions[-1], bins=200, range=(-1.75, 1.75), density=True)
    # x = np.linspace(-1.75, 1.75, 1000)
    # plt.plot(x, levy_stable.pdf(x, alpha=alpha, beta=0, loc=mu, scale=sigma), '--', lw=2, label="L\'evy")
    # from scipy.stats import norm
    # plt.plot(x, norm.pdf(x, loc=0, scale=np.std(states.emissions[-1])), '--', lw=2, label='Gaussian')
    # plt.legend()
    # plt.tight_layout()
    # plt.show()
    # exit()
    #############################################################

    # states.transition_autocorrelation(plot=True)
    # states.calculate_hurst()
    # # print(states.fit_params)
    # for i in states.hurst:
    #     print(i.mean())
    # exit()
    print(states.transition_matrix)
    exit()
    import itertools

    dwells = [[] for _ in range(states.nstates)]
    for solute in range(states.nsolute):
        for state, group in itertools.groupby(states.state_sequence[:, solute]):
            dwells[state].append(len(list(group)))

    all_dwells = []
    for i in dwells:
        all_dwells += i

    plt.hist(all_dwells, bins=22, range=(3, 25))
    plt.show()
    exit()
    #
    # for i in range(8):
    #     plt.hist(dwells[i], bins=50)
    #     plt.show()
    # exit()
    #
    # emissions = []
    # for i in range(states.nstates):
    #     for j in range(states.nstates):
    #         if i != j:
    #             emissions += states.emissions[i][j]
    #
    # plt.hist(emissions, bins=50, range=(-.75, .75))
    # plt.show()
    # exit()

    chains = Chain(states.count_matrix, states.fit_params, hurst_parameters=states.hurst,
                   emission_function=levy_stable)
    # chains.generate_realizations(24, 10000, bound=states.maximum_emission())
    chains.generate_realizations(24, 10000, bound=1.0)
    chains.calculate_msd()
    chains.plot_msd()
