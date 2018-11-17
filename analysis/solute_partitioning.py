#! /usr/bin/env python

import argparse
import mdtraj as md
from LLC_Membranes.llclib import physical, topology
from LLC_Membranes.analysis import Atom_props
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from scipy import spatial
from pymbar import timeseries
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
    parser.add_argument('-spline', action="store_true", help='Trace center of pore with a spline instead of'
                                                             'approximating each center as a single point.')
    parser.add_argument('-pr', '--pore_radius', default=0.5, type=float, help='Max distance from pore center a water '
                        'molecule can exist in order to be counted as inside the pore')
    parser.add_argument('-b', '--buffer', default=0, type=float, help='z-distance (nm) to cut off from top and bottom '
                                                                      'of unit cell when calculating COM locations')

    parser.add_argument('--load', action="store_true")
    parser.add_argument('--savename', default='solute_partitioning.pl')

    parser.add_argument('-boot', '--nboot', default=200, type=int, help='Number of bootstrap trials')
    parser.add_argument('--single_frame', action='store_true', help='Specify this flag in order to analyze a single'
                                                                    '.gro file. No statistics will be generated')
    parser.add_argument('--tcl', action='store_true', help='Create .tcl file that create representation of water within'
                                                           'pore region. Use vmd *.gro -e *.tcl. Only designed for use'
                                                           'with a single frame')

    return parser


def make_groups(t, pos, pore_atoms, tail_atoms, bounds, natoms, mpl):
    """
    :param t: mdtraj trajectory object
    :param pore_atoms: atoms that define the pore region, list
    :param tail_atoms: atoms that define the tail region, list
    :param bounds: only include atoms within these z boundaries. list:[lower_bound, upper_bound]
    :param natoms: number of atoms in monomer residue
    :param mpl: number of monomers per layer
    :return:
    """

    all_pore_atoms = []
    all_tail_atoms = []
    water = []  # oxygen atom of water molecule
    other = []  # all other atoms within bounds. They will be needed for a weight percent calculation
    mass = 0
    for a in t.topology.atoms:
        if bounds[0] <= pos[a.index, 2] <= bounds[1]:
            mass += Atom_props.mass[a.name]
            if a.residue.name == 'HOH' and a.name == 'O':
                water.append(a.index)
            elif a.name in pore_atoms:
                all_pore_atoms.append(a.index)
            elif a.name in tail_atoms:
                all_tail_atoms.append(a.index)

    # all_pore_atoms will be used to make a spline which traces through the pores. To do that properly, only full layers
    # can be included in the list.

    # layer_number = [i // (mpl*natoms) for i in all_pore_atoms]
    #
    # pore_atoms_keep = []
    # for i in range(len(all_pore_atoms)):
    #     if layer_number.count(layer_number[i]) == int(mpl*len(pore_atoms)):
    #         pore_atoms_keep.append(all_pore_atoms[i])
    #     else:
    #         other.append(all_pore_atoms[i])

    return all_pore_atoms, all_tail_atoms, water, mass


def restrict_spline(full_spline, bounds):
    """
    Restrict a pore-tracing spline to be within two z boundaries
    :param full_spline: A full spline
    :param bounds: z boundaries
    :return: restricted spline to within z-boundaries
    """

    restricted = []
    for i in range(full_spline.shape[0]):
        if bounds[0] <= full_spline[i, 2] <= bounds[1]:
            restricted.append(full_spline[i, :])

    return np.array(restricted)


def water_content(pos, ref_pos, r, write=False):
    """
    Find number of water molecules in region
    :param pos: positions of all water molecules (excluding those outside boundaries defined by args.bounds)
    :param ref_pos: positions defining the region
    :param r: distance from reference positions a point can be for it to be included as a part of the region
    :return: number of water molecules in region
    """

    n = 0
    ndx = []
    tree = spatial.cKDTree(ref_pos)
    for i in range(pos.shape[0]):
        d = tree.query(pos[i, :])[0]
        if d < r:
            n += 1
            ndx.append(i + 1)

    if write:
        return n, ndx
    else:
        return n


def bootstrap(data, tau, nboot):
    """
    :param data: equilibrated data from a timeseries
    :param tau: autocorrelation time (number of points in timeseries between uncorrelated samples)
    :param nboot: number of bootstrap trials
    :return: average and standard deviation
    """

    tau = int(tau)
    nsub = data.shape[0] // tau  # number of uncorrelated subtrajectories to break data into

    # divide data into independent subtrajectories
    subtrajectories = np.zeros([nsub, tau])
    for i in range(nsub):
        if i == 0:
            subtrajectories[i, :] = data[-tau:]
        else:
            subtrajectories[i, :] = data[-(i+1)*tau:-i*tau]

    choices = np.linspace(0, tau - 1, tau, dtype=int)  # indices of each subtrajectory
    boot = np.zeros([nboot])
    for b in range(nboot):
        trial = np.zeros([nsub])  # the statistic will be the average of each independent subtrajectory
        for s in range(nsub):
            ndx = np.random.choice(choices, size=tau, replace=True)
            trial[s] = np.mean(subtrajectories[s, choices])
        boot[b] = np.mean(trial)

    return np.mean(boot), np.std(boot)


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

        self.ids = [a.name for a in self.t.topology.atoms]
        self.res = [a.residue.name for a in self.t.topology.atoms]

        self.pore_atoms = pore_atoms
        self.npores = npores
        self.pore_centers = None
        self.pore_water = []
        self.tail_water = []

        if residue == 'SOL':
            for a in self.t.topology.atoms:
                if a.residue.name == 'HOH':  # workaround for mdtraj
                    if a.name == 'O':
                        a.name = 'OW'
                    elif a.name == 'H1':
                        a.name = 'HW1'
                    elif a.name == 'H2':
                        a.name = 'HW2'
            for a in self.t.topology.atoms:
                if a.residue.name == 'HOH':
                    a.residue.name = 'SOL'

        self.residue = topology.Residue(residue)
        self.residue_indices = np.array([a.index for a in self.t.topology.atoms if a.residue.name == residue])
        residue_atom_names = [a.name for a in self.t.topology.atoms if a.residue.name == residue]
        masses = [self.residue.mass[x] for x in residue_atom_names[:self.residue.natoms]]
        self.com = physical.center_of_mass(self.pos[:, self.residue_indices, :], masses)

    def locate_pore_centers(self, spline=False):

        # find pore centers
        pore_atoms = [a.index for a in self.t.topology.atoms if a.name in self.pore_atoms]
        if spline:
            print('Creating pore splines')
            self.pore_centers = physical.trace_pores(self.pos[:, pore_atoms, :], self.t.unitcell_vectors, 20)
        else:
            self.pore_centers = physical.avg_pore_loc(self.npores, self.pos[:, pore_atoms, :])

    def partition(self, r, buffer=0):
        """ Partition solute residue into tail and pore region

        :param r: pore radius, outside of which atoms will be considered in the tail region
        :param buffer: z distance (nm) to cut out from top and bottom of membrane (in cases where there is a water gap)

        :type r: float
        :type buffer: float

        """

        # to help determine buffer
        # plt.hist(self.com[0, :, 2], bins=50)
        # plt.show()

        for i in tqdm.tqdm(range(self.t.n_frames)):
            pore = []

            if buffer > 0:
                xy_positions = self.com[i, (self.com[i, :, 2] > buffer) & (self.com[i, :, 2] <
                                        self.t.unitcell_vectors[i, 2, 2] - buffer), :2]
            else:
                xy_positions = self.com[i, :, :2]

            for p in range(self.npores):
                d = np.linalg.norm(xy_positions - self.pore_centers[i, p, :], axis=1)
                pore += np.where(d <= r)[0].tolist()

            self.pore_water.append(np.unique(pore).tolist())

            tails = np.full(xy_positions.shape[0], True, dtype=bool)  # array of True. Booleans are fast
            tails[self.pore_water[i]] = False  # set pore indices to False
            self.tail_water.append(np.where(tails)[0])  # list of tail indices where tails = True

            # vvv slow vvv
            #self.tail_water.append([x for x in np.arange(0, self.com.shape) if x not in self.pore_water[i]])

    def plot(self, resname='Water'):

        ntail_water = [len(x) for x in self.tail_water]
        npore_water = [len(x) for x in self.pore_water]

        plt.plot(self.t.time/1000, ntail_water, color='xkcd:blue', label='Tail %s' % resname)
        plt.plot(self.t.time/1000, npore_water, color='xkcd:orange', label='Pore %s' % resname)
        plt.xlabel('Time (ns)')
        plt.ylabel('Number of water molecules')
        plt.legend(fontsize=14)

        plt.tight_layout()

        plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if not args.load:
        # heavy calcuations
        sys = System(args.gro, args.pore_atoms, args.residue, traj=args.traj, begin=args.begin, end=args.end,
                     skip=args.skip)

        sys.locate_pore_centers(spline=args.spline)

        with open(args.savename, "wb") as f:
            pickle.dump(sys, f)
    else:

        with open(args.savename, "rb") as f:
            sys = pickle.load(f)

    print('Calculating solute partition by frame')
    sys.partition(args.pore_radius, buffer=args.buffer)
    sys.plot(resname='Water')

