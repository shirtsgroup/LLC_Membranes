#! /usr/bin/env python

"""
Calculate H-bonds in system.

Some definitions:
    donor (D): atom covalently bonded to hydrogen
    acceptor (A) : atom which 'accepts' hydrogen bond

A diagram of an hbond:

    D--H - - A

Criterion:
    Distance between D and A below some distance
    Angle between DHA less than some cut-off

Exact values for criterion are left up to the user.
"""


import argparse
import mdtraj as md
import numpy as np
import os
import tqdm
import matplotlib.pyplot as plt
import pickle
from LLC_Membranes.llclib import file_rw, topology, physical


def initialize():

    parser = argparse.ArgumentParser(description='Use geometric criteria to identify hydrogen bonds')

    parser.add_argument('-t', '--traj', default='wiggle.trr', type=str, help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='End frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Start frame')
    parser.add_argument('-sk', '--skip', default=1, type=int, help='Skip every skip frames')
    parser.add_argument('-r', '--residues', nargs='+', default=['HII'], help='Residues to include in h-bond search')
    parser.add_argument('-a', '--atoms', action='append', nargs='+', help='Atoms to include for each '
                        'residue. Each list of atoms must be passed with a separate -a flag for each residue in'
                        'args.residues')
    parser.add_argument('-d', '--distance', default=.35, type=float, help='Maximum distance between acceptor and donor'
                                                                         ' atoms')
    parser.add_argument('-angle', '--angle_cut', default=30, type=float, help='Maximum DHA angle to be considered an H-bond')
    parser.add_argument('-l', '--load', default=False, help='Load pickled system object')
    parser.add_argument('-acc', '--acceptors', default=False, nargs='+', help='If you only want hbonds with acceptor'
                        'atoms of a residue, use this flag. 1 for this condition, otherwise 0. Input as a list'
                        ' of 0 and 1 in the same order as args.residues')
    parser.add_argument('-donors', '--donors', default=False, nargs='+', help='If you only want hbonds with donor'
                        'atoms of a residue, use this flag. 1 for this condition, otherwise 0. Input as a list'
                        ' of 0 and 1 in the same order as args.residues')
    parser.add_argument('-s', '--savename', default='hbonds.pl', type=str, help='Name of pickled object to save after'
                                                                                'calcualtions are complete')
    parser.add_argument('-noshow', '--noshow', action="store_true", help='Do not display plots')

    return parser


class System(object):

    def __init__(self, traj, gro, begin=0, end=-1, skip=1, t=None):
        """ Load in the trajectory

        :param traj: GROMACS trajectory file (.xtc or .trr)
        :param gro: GROMACS coordinate file (.gro)
        :param begin: First frame to analyze
        :param end: Last frame to analyze
        :param skip: Only analyze every `skip` frames
        :param t: mdtraj trajectory object. If this is provided, traj and gro will not be loaded

        :type traj: str
        :type gro: str
        :type begin: int
        :type end: int
        :type skip: int
        :type t: object
        """

        # print('Generating bond list...', end="", flush=True)
        # self.top = Topology(top, xlink=xlink, xlink_topology=xlink_topology, xlink_residue=xlink_residue)
        # print('Done!')

        if t is None:
            print('Loading trajectory...', end="", flush=True)
            self.t = md.load(traj, top=gro)[begin:end:skip]
            print('Done!')
        else:
            self.t = t

        self.pos = self.t.xyz  # positions of all atoms
        self.hbonds = []  # will hold h-bonds for each frame [D, H, A, angle]
        self.dt = self.t.time[1] - self.t.time[0]
        self.box = self.t.unitcell_vectors

        # for hbond_pairing
        # self.nwater = 0
        # self.nres = 0
        # self.water_numbers = None
        # self.residue_numbers = None

        # system definitions
        self.names = []
        for a in self.t.topology.atoms:
            # workaround for mdtraj since standard_names is not an option with .gro topologies.
            if a.residue.name == 'HOH':
                if a.name == 'O':
                    self.names.append('OW')
                elif a.name == 'H1':
                    self.names.append('HW1')
                elif a.name == 'H2':
                    self.names.append('HW2')
            else:
                self.names.append(a.name)

        names = topology.fix_names(gro)  # rename atoms because mdtraj screws it up in some cases.
        for i, a in enumerate(self.t.topology.atoms):
            a.name = names[i]

        self.residues = [a.residue.name for a in self.t.topology.atoms]

        # definitions
        self.donor_atoms = ['H']
        self.acceptor_atoms = ['O', 'N']  # atoms that can be acceptors

        # initialize lists
        self.D = []
        self.H = []
        self.A = []

        self.donor_acceptor_matrix = None
        self.atom_to_matrix_index = {}
        self.matrix_to_atom_index = {}

        # Now, specify water with the -r flag
        # if not exclude_water:
        #
        #     # all H's are potential donors
        #     self.H = [a.index for a in self.t.topology.atoms if a.residue.name == 'HOH' and a.element.symbol == 'H']
        #     # get the index of the atoms bonded to each H (all oxygens)
        #     self.D = [self.top.bonds[x][0] for x in self.H]  # assumes only one bond to H as it should
        #     # all oxygens are also potential acceptors
        #     self.A = [a.index for a in self.t.topology.atoms if a.residue.name == 'HOH' and a.element.symbol == 'O']

    def set_eligible(self, res, atoms, acceptor_only=False, donors_only=False):
        """ Define atoms that will be included in h-bond calculation. Make sure `res` is updated with proper
        annotations or no atoms will be selected!

        :param res: residue to include in calculation
        :param atoms: atoms from residue to include in calculation. Use 'all' to include all potential hbonding atoms
        :param acceptor_only: restrict residue so it can only accept hydrogen bonds
        :param donor_only: restrict residue so it can only donate hydrogen bonds

        :type res: str
        :type atoms: str or list
        :type acceptor_only: bool
        :type donors_only: bool
        """

        residue = topology.Residue(res)

        restrict = False
        if type(atoms) is list and atoms[0] != 'all':
            restrict = True
        elif type(atoms) is not list and atoms != 'all':
            restrict = True

        if restrict:
            residue.hbond_H = [h for h in residue.hbond_H if h in atoms]
            residue.hbond_D = [d for d in residue.hbond_D if d in atoms]
            residue.hbond_A = [a for a in residue.hbond_A if a in atoms]

        if acceptor_only:  # if we only want this residue to recieve hbonds, don't include its donor atom or H
            residue.hbond_H = []
            residue.hbond_D = []

        if donors_only:
            residue.hbond_A = []

        for a in self.t.topology.atoms:
            if a.residue.name == res:
                if a.name in residue.hbond_H:
                    self.H.append(a.index)
                if a.name in residue.hbond_D:
                    self.D.append(a.index)
                if a.name in residue.hbond_A:
                    self.A.append(a.index)

    def identify_hbonds(self, cut, angle):
        """ Identify hydrogen bonds based on geometric criteria. If the angle formed by D-H--A and the distance between
        donor and acceptor are less than the cut-offs, then we consider a hydrogen bond to exist

        :param cut: cut-off distance (nm)
        :param angle: cut-off angle (degrees)

        :type cut: float
        :type angle: float
        """

        # narrow list by doing a distance search

        Dlen = len(self.D)
        Alen = len(self.A)
        d = np.zeros([self.pos.shape[0], Dlen * Alen])

        # distance search can be sped up with cKDtree
        if Dlen >= Alen:
            L = Dlen
            print('Calculating distances...')
            for t in tqdm.tqdm(range(d.shape[0]), unit='frames'):
                for i in range(len(self.A)):
                    minimum_image_distance = physical.minimum_image_distance(self.pos[t, self.D, :] -
                                                                             self.pos[t, self.A[i], np.newaxis, :],
                                                                             self.box[t, ...])
                    d[t, i * Dlen:(i + 1) * Dlen] = np.linalg.norm(minimum_image_distance, axis=1)
        else:
            L = Alen
            print('Calculating distances...')
            for t in tqdm.tqdm(range(d.shape[0]), unit='frames'):
                for i in range(len(self.D)):
                    minimum_image_distance = physical.minimum_image_distance(self.pos[t, self.A, :] -
                                                                             self.pos[t, self.D[i], np.newaxis, :],
                                                                             self.box[t, ...])
                    d[t, i * Alen:(i + 1) * Alen] = np.linalg.norm(minimum_image_distance, axis=1)

        print('Narrowing eligible atoms and calculating angles')

        for i in tqdm.tqdm(range(d.shape[0])):

            indices = np.where(d[i, :] < cut)[0]  # indices where distance is below cutoff

            distance_eligible = indices[np.nonzero(d[i, indices])]  # narrow to only nonzero values

            if Dlen >= Alen:
                Aindex = np.array(self.A)[distance_eligible // L]  # indices of distance eligible acceptor atoms
                Dindex = np.array(self.D)[distance_eligible % L]  # indices of distance eligible donor atoms
                Hindex = np.array(self.H)[distance_eligible % L]  # H atoms attached to eligible donors
            else:
                Aindex = np.array(self.A)[distance_eligible % L]  # indices of distance eligible acceptor atoms
                Dindex = np.array(self.D)[distance_eligible // L]  # indices of distance eligible donor atoms
                Hindex = np.array(self.H)[distance_eligible // L]  # H atoms attached to eligible donors

            self.hbonds.append(np.reshape(np.concatenate((Dindex, Hindex, Aindex)), (3, len(Aindex))))

            if self.hbonds[i].size > 0:  # don't calculate vectors if there aren't any hbonds this frame
                # calculate vectors
                v = np.zeros([2, len(Aindex), 3])

                v[0, ...] = physical.minimum_image_distance(self.pos[i, self.hbonds[i][1, :], :] -
                                                            self.pos[i, self.hbonds[i][0, :], :], self.box[i, ...])  # D-H vectors
                v[0, ...] /= np.linalg.norm(v[0, ...], axis=1)[:, np.newaxis]  # normalize (need to, to get correct angle)
                v[1, :] = physical.minimum_image_distance(self.pos[i, self.hbonds[i][2, :], :] -
                                                          self.pos[i, self.hbonds[i][1, :], :], self.box[i, ...])  # H-A vectors
                v[1, ...] /= np.linalg.norm(v[1, ...], axis=1)[:, np.newaxis]  # normalize

                # calculate angles

                dot = np.zeros([v.shape[1]])
                for j in range(dot.shape[0]):
                    dot[j] = np.dot(v[0, j, :], v[1, j, :])

                a = np.arccos(dot) * (180 / np.pi)  # convert to degrees

                self.hbonds[i] = np.delete(self.hbonds[i], np.where(a > angle)[0], axis=1)
                self.hbonds[i] = np.concatenate((self.hbonds[i], a[a < angle][np.newaxis, :]), 0)

    def plot_hbonds(self, show=True, save=True, savename='hbonds.png'):
        """ Plot the total number of hydrogen bonds per frame as a function of time

        :param show: display plot
        :param save: save plot under savename
        :param savename: name under which to save plot

        :type show: bool
        :type save: bool
        :type savename: str
        """

        n = [a.shape[1] for a in self.hbonds]

        # ethers = ['O', 'O1', 'O2']
        # carb = ['O3', 'O4']
        # nether = np.zeros([self.t.n_frames])
        # ncarb = np.zeros([self.t.n_frames])
        #
        # for i, x in enumerate(self.hbonds):
        #
        #     if x.shape[0] == 4:
        #
        #         acceptors = [self.names[int(x[2, j])] for j in range(x.shape[1])]
        #
        #         for k in acceptors:
        #             if k in ethers:
        #                 nether[i] += 1
        #             elif k in carb:
        #                 ncarb[i] += 1
        #
        # print(sum(n))
        # print(np.sum(ncarb) / np.sum(nether))

        # plt.plot(self.t.time / 1000, nether / ncarb)

        print("Average Hydrogen Bonds Per Frame: % .2f" % np.mean(n))
        plt.plot(self.t.time / 1000, n)
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Number of hydrogen bonds', fontsize=14)
        plt.tight_layout()
        if save:
            plt.savefig(savename)
        if show:
            plt.show()

    def number_water_molecules(self):
        """
        A specific form of number_residues function (below)
        :return:
        """
        water_numbers = {}
        nwater = 0

        atoms = 1
        for a in self.t.topology.atoms:
            if a.residue.name == 'HOH':
                water_numbers[a.index] = nwater
                if atoms % 3 == 0:
                    nwater += 1
                atoms += 1

        return water_numbers

    def number_residues(self, res):

        residue_numbers = {}
        nres = 0

        residue = Residue(res)
        atoms = 1
        for a in self.t.topology.atoms:
            if a.residue.name == res:
                residue_numbers[a.index] = nres
                if atoms % residue.natoms == 0:
                    nres += 1
                atoms += 1

        return residue_numbers, nres

    def save_hbonds(self, name='hbonds.pl'):

        with open(name, 'wb') as f:
            pickle.dump(self.hbonds, f)

    def hbond_matrix(self):

        unique_atoms = []
        for i in range(len(self.hbonds)):
            unique_atoms.append(np.unique(self.hbonds[i][::2, :]))

        unique = np.unique(np.concatenate(unique_atoms))

        for i in range(unique.size):
            self.atom_to_matrix_index[int(unique[i])] = i
            self.matrix_to_atom_index[i] = int(unique[i])

        self.donor_acceptor_matrix = np.zeros([len(self.hbonds), unique.size])

        for t in range(len(self.hbonds)):
            hbond_frame = self.hbonds[t]
            for i in range(hbond_frame.shape[1]):
                self.donor_acceptor_matrix[t, self.atom_to_matrix_index[hbond_frame[0, i]]] = hbond_frame[2, i]


class Residue(object):

    def __init__(self, name, xlink=False, xlink_topology='assembly.itp'):

        self.name = name

        script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

        self.is_ion = False
        # check if residue is an ion
        with open('%s/../top/topologies/ions.txt' % script_location) as f:
            ions = []
            for line in f:
                if line[0] != '#':
                    ions.append(str.strip(line))

        if name in ions:
            self.is_ion = True
            self.natoms = 1

        else:

            try:
                if xlink:
                    f = open('%s' % xlink_topology, 'r')
                else:
                    f = open('%s.itp' % name, 'r')
            except FileNotFoundError:
                try:
                    f = open('%s/../top/topologies/%s.itp' % (script_location, name), 'r')
                except FileNotFoundError:
                    sys.exit('No topology %s.itp found' % name)

            itp = []
            for line in f:
                itp.append(line)

            f.close()

            atoms_index = 0
            while itp[atoms_index].count('[ atoms ]') == 0:
                atoms_index += 1

            self.indices = {}  # key = atom name , value = index
            self.names = {}  # key = index, value = atom name

            atoms_index += 2
            i = 0
            while itp[i + atoms_index] != '\n':
                data = itp[i + atoms_index].split()
                self.indices[data[4]] = int(data[0])
                self.names[int(data[0])] = data[4]
                i += 1

            self.natoms = len(self.names)

            # find the bonds section
            bonds_index = 0
            while itp[bonds_index].count('[ bonds ]') == 0:
                bonds_index += 1
            bonds_index += 2

            bonds = []
            while itp[bonds_index] != '\n':
                bond_data = str.split(itp[bonds_index])[:2]
                bonds.append([int(bond_data[0]), int(bond_data[1])])
                bonds_index += 1

            bonds = np.array(bonds)  # converted to numpy array for use with np.where

            self.bonds = {}
            for i in range(self.natoms):
                self.bonds[i] = []

                # involvement = [list(x) for x in bonds if i + 1 in x]  # old (slow) way of doing this

                x = np.where(bonds == i + 1)
                involvement = [bonds[x[0][i]] for i in range(len(x[0]))]

                for pair in involvement:
                    atom = [x - 1 for x in pair if x != (i + 1)][0]
                    self.bonds[i].append(atom)


if __name__ == "__main__":

    args = initialize().parse_args()

    if not args.acceptors:
        args.acceptors = [False for i in args.residues]
    else:
        args.acceptors = [bool(int(i)) for i in args.acceptors]

    if not args.donors:
        args.donors = [False for i in args.residues]
    else:
        args.donors = [bool(int(i)) for i in args.donors]

    if args.load:

        sys = file_rw.load_object(args.load)

    else:
        # workaround for argparse. If default value is set, it is always included in the list with action='append'
        if not args.atoms:
            args.atoms = [['all'] for r in args.residues]  # a default value

        while len(args.atoms) != len(args.residues):
            args.atoms.append(['all'])

        sys = System(args.traj, args.gro, begin=args.begin, end=args.end, skip=args.skip)

        for i, r in enumerate(args.residues):
            sys.set_eligible(r, args.atoms[i], acceptor_only=args.acceptors[i], donors_only=args.donors[i])

        sys.identify_hbonds(args.distance, args.angle_cut)
        file_rw.save_object(sys, '%s' % args.savename)

    if not args.noshow:
        sys.plot_hbonds()
