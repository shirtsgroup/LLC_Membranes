#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
from LLC_Membranes.llclib import physical, topology, transform
import sys
import tqdm
from scipy.sparse import lil_matrix
import pickle
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

    # atom selection
    parser.add_argument('-r', '--residue', default=None, help='Residue to calculate coordination number with respect '
                                                               'to')
    parser.add_argument('-rc', '--coordinated_residue', default=None, help='Name of residue coordinated to residue_')
    parser.add_argument('-ac', '--coordinated_atoms', default=None, nargs='+', help='Name of residue coordinate to '
                        'residue')
    parser.add_argument('-a', '--atoms', default=None, nargs='+', help='Name of atoms to calculate correlation '
                        'function with respect to. The center of mass will be used')
    parser.add_argument('-ta', '--type', default=None, help='Element name of atoms that you want coordination number of')
    parser.add_argument('-tc', '--coordinated_type', default=None, help='Element name of coordinated atoms')

    # coordination calculation parameters
    parser.add_argument('-cut', default=0.25, type=float, help='Maximum distance between pairs where they are considered'
                                                              'coordinated (nm)')

    # saving options
    parser.add_argument('-s', '--savename', default='coordination.pl', help='Name under which to save distance '
                                                                                 'array')
    parser.add_argument('-l', '--load', action="store_true", help='Load distance array saved with name '
                                                                             'passed to save_distances option.')

    parser.add_argument('-bins', nargs='+', default=100, type=int, help='Integer or array of bin values. If more than'
                        'one value is used, order the inputs according to the order given in args.axis')
    parser.add_argument('-nocom', action='store_true', help='Calculate coordination based on pairwise distance between'
                                                            'all atoms, even those that are a part of the same residue.')

    return parser


class System(object):

    def __init__(self, traj, gro, residue=None, coordinated_residue=None, atoms=None, coordinated_atoms=None, type=None,
                 ctype=None, begin=0, end=-1, skip=1, com=True):
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
        """

        print("Loading trajectory...", end='', flush=True)
        self.t = md.load(traj, top=gro)[begin:end:skip]
        print("Done!")

        names = topology.fix_names(gro)  # rename atoms because mdtraj screws it up in some cases.
        for i, a in enumerate(self.t.topology.atoms):
            a.name = names[i]

        self.time = self.t.time / 1000  # time in nanoseconds

        # get locations of atoms of interest and coordinating atoms
        # calculate centers of mass for certain groups. Decision-making in this regard is made in self.narrow_atoms()

        if atoms is not None:
            if atoms[0] == 'all':
                atoms = self.all_atoms(res=residue)

        if coordinated_atoms is not None:
            if coordinated_atoms[0] == 'all':
                coordinated_atoms = self.all_atoms(res=coordinated_residue)

        self.com, self.com_map = self.narrow_atoms(atoms, residue, type, com=com)
        self.com_coordinated, self.com_coordinated_map = self.narrow_atoms(coordinated_atoms, coordinated_residue,
                                                                           ctype, coordination=True, com=com)

        # relate indices to
        self.names = [a.name for a in self.t.topology.atoms]
        self.residues = [a.residue.name for a in self.t.topology.atoms]

        self.distances = None
        self.ncoord = None

    def all_atoms(self, res=None):

        if res is not None:

            r = topology.Residue(res)

            return [x for x in r.names.values()]

    def narrow_atoms(self, atoms, residue, type, coordination=False, com=True):

        if atoms is not None or type is not None:

            if residue is not None:

                if type is not None:
                    # get names of all atoms with the appropriate type
                    atoms = set([a.name for a in self.t.topology.atoms if a.element.symbol == type])

                # get indices of all atoms in the system that make up the atoms list
                atom_indices = [a.index for a in self.t.topology.atoms if a.residue.name == residue and
                                a.name in atoms]

                if type is not None:
                    # if a residue is specified with an atom type, assume you want the locations of each atom of that
                    # type within each residue.

                    return self.t.xyz[:, atom_indices, :], topology.map_atoms(atom_indices)

                else:
                    # If a residue is specified with atoms, assume that the user wants to isolate calculations to the
                    # center of mass of a group of atoms within a particular residue

                    if com:

                        residue = topology.Residue(residue)
                        atom_mass = [residue.mass[v] for v in residue.mass.keys() if
                                     v in atoms]  # mass of atoms of interest

                        if len(atom_mass) > 1:

                            return physical.center_of_mass(self.t.xyz[:, atom_indices, :], atom_mass), \
                                   topology.map_atoms(atom_indices, len(atom_mass))

                        else:

                            self.t.xyz[:, atom_indices, :], topology.map_atoms(atom_indices)

                    else:

                        return self.t.xyz[:, atom_indices, :], topology.map_atoms(atom_indices)

            else:

                if type is not None:
                    # get indices of any atoms whose element = type
                    atom_indices = [a.index for a in self.t.topology.atoms if a.element.symbol == type]

                else:
                    # get indices of any atoms in the "atoms" list
                    atom_indices = [a.index for a in self.t.topology.atoms if a.name in atoms]

            return self.t.xyz[:, atom_indices, :], topology.map_atoms(atom_indices)

        else:

            if residue is not None:
                # calculate the center of mass of the residue based on all atoms
                # First get indices of all atoms in the residue
                atom_indices = [a.index for a in self.t.topology.atoms if a.residue.name == residue]
                res = topology.Residue(residue)
                atom_mass = [v for v in res.mass.values()]  # mass of each atom in an individual residue

                return physical.center_of_mass(self.t.xyz[:, atom_indices, :], atom_mass), \
                       topology.map_atoms(atom_indices, len(atom_mass))

            else:
                # if you forget a flag, exit the program with a descriptive error
                if coordination:
                    sys.exit('You must supply at least a residue (-cr / --coordinated_residue) or an atom '
                             '(-ca / --coordinated_atoms)')
                else:
                    sys.exit('You must supply at least a residue (-r / --residue) or an atom name (-a / --atoms)')

    def distance_search(self, cut=0.31):
        """ Find all minimum image pairwise distances

        :param cut: maximum distance

        :type cut: float

        :return:
        """

        # initialize array to hold all pairwise distances
        self.distances = np.zeros([self.t.n_frames], dtype=object)

        print('Calculating minimum image distances!')
        for t in tqdm.tqdm(range(self.distances.shape[0])):
            self.distances[t] = lil_matrix((self.com.shape[1], self.com_coordinated.shape[1]))  # sparse matrices are 2D
            for i in range(self.com.shape[1]):
                xyz_distances = self.com_coordinated[t, ...] - self.com[t, i, :]
                min_distances = physical.minimum_image_distance(xyz_distances, self.t.unitcell_vectors[t, ...])
                euclidean_dist = np.linalg.norm(min_distances, axis=1)
                under_cut = np.where(euclidean_dist < cut)[0]
                self.distances[t][i, under_cut] = euclidean_dist[under_cut]

    def n_coordinated(self, res=None, atom_groups=None, plot=True):
        """ Plot number of atoms coordinated to residue vs. time

        :param res: a list of residue names. This function will plot coordinated atoms belonging to specific residues
        separately.
        :param atom_groups: a list of lists of atom names. Can be used to distinguish between coordination with
        different groups of atoms within the same residue. If atom names are the same between residues, also pass
        the appropriate residue name for each group of atoms to res

        :type res: list
        :type atom_groups: list

        :return:
        """

        if atom_groups is not None:

            self.ncoord = np.zeros([len(atom_groups), self.t.n_frames, self.com.shape[1]])

            if res is not None:

                print('Organizing atoms based on atom type and residue name...')
                for t in tqdm.tqdm(range(self.t.n_frames)):
                    for i in range(self.ncoord.shape[2]):
                        atom_types = [self.names[self.com_coordinated_map[j][0]] for j in
                                      np.nonzero(self.distances[t][i, :])[1]]  # residue name
                        residues = [self.residues[self.com_coordinated_map[j][0]] for j in
                                    np.nonzero(self.distances[t][i, :])[1]]

                        grp = []
                        for ndx, k in enumerate(atom_types):
                            atoms = np.array([a.count(k) for a in atom_groups])
                            resi = np.array([a.count(residues[ndx]) for a in res])
                            try:
                                grp.append(np.where(atoms + resi == 2)[0][0])
                            except IndexError:
                                pass

                        for g in grp:
                            self.ncoord[g, t, i] += 1

            else:

                print('Organizing atoms based on atom type...')
                for t in tqdm.tqdm(range(self.t.n_frames)):
                    for i in range(self.ncoord.shape[2]):
                        atom_types = [self.names[self.com_coordinated_map[j][0]] for j in
                                      np.nonzero(self.distances[t][i, :])[1]]  # residue name

                        grp = []
                        for ndx, k in enumerate(atom_types):
                            atoms = [a.count(k) for a in atom_groups]

                            try:
                                grp.append(atoms.index(1))
                            except IndexError:
                                pass

                        for g in grp:
                            self.ncoord[g, t, i] += 1

            # for i in range(ncoord.shape[0]):
            #
            #     plt.plot(self.time, ncoord[i, ...].mean(axis=1), label='%s of %s' % (atom_groups[i], res[i]),
            #              linewidth=2)
            colors = ['blue', 'green', 'orange']
            for j in np.random.randint(self.ncoord.shape[2], size=1):
                for i in range(self.ncoord.shape[0]):

                    plt.plot(self.time, self.ncoord[i, :, j], label='%s of %s' % (atom_groups[i], res[i]),
                             linewidth=2, color=colors[i])

            #plt.legend()

        else:
            self.ncoord = np.zeros([self.t.n_frames, self.com.shape[1]])

            for i in tqdm.tqdm(range(self.ncoord.shape[1])):
                self.ncoord[:, i] = [len(np.nonzero(self.distances[t][i, :])[1]) for t in range(self.t.n_frames)]
                #self.ncoord[:, i] = [(self.distances[t] != 0).sum(1).mean() for t in range(self.t.n_frames)]

            # for i in range(self.com.shape[1]):
            #     print(i)
            #     plt.plot(self.ncoord[:, i])
            #     plt.show()
            #
            # exit()
            print('Average coordinated molecules per frame: %.2f' % self.ncoord.mean())
            #plt.plot(self.time, self.ncoord.mean(axis=1))
            plt.plot(self.time, self.ncoord.sum(axis=1))

        if plot:

            plt.xlabel('Time (ns)', fontsize=14)
            plt.ylabel('Number of coordinated molecules', fontsize=14)
            plt.gcf().get_axes()[0].tick_params(labelsize=14)
            plt.tight_layout()

            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if not args.load:

        if args.nocom:
            com = False
        else:
            com = True

        system = System(args.traj, args.gro, atoms=args.atoms, coordinated_atoms=args.coordinated_atoms,
                        residue=args.residue, coordinated_residue=args.coordinated_residue, type=args.type,
                        ctype=args.coordinated_type, begin=args.begin, end=args.end, skip=args.skip, com=com)

        system.distance_search(cut=args.cut)  # calculate pairwise distance between all points in self.com and self.com_coordinated

        system.n_coordinated(plot=True)

        print((system.ncoord.flatten() > 0).sum() / (system.ncoord.shape[0] * system.ncoord.shape[1]))

        # with open(args.savename, 'wb') as f:
        #     pickle.dump(system, f)

    else:
        print('Loading pickled object!...', end='', flush=True)
        system = pickle.load(open(args.savename, "rb"))
        print('Done!')

        plt.plot(system.time, system.ncoord.sum(axis=1))
        plt.show()
        exit()
        print((system.ncoord.flatten() > 0).sum() / (system.ncoord.shape[0] * system.ncoord.shape[1]))

        # system.n_coordinated(plot=True)

    #system.plot(res=['HII', 'HII', 'HOH'], atom_groups=[['O3', 'O4'], ['O', 'O1', 'O2'], ['O']])
    # system.n_coordinated(plot=True)
    #
    # for i in np.nonzero(system.distances[-1][0, :])[1]:
    #     print('%s -- %s' % (system.com_map[0], system.com_coordinated_map[i]))
