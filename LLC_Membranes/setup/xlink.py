#!/usr/bin/env python

import mdtraj as md
import argparse
import time
import os
import numpy as np
from LLC_Membranes.setup.xlink_schemes import XlinkReaction
from LLC_Membranes.setup.add_dummies import add_dummies
from LLC_Membranes.setup.gentop import SystemTopology
from LLC_Membranes.setup.genmdp import SimulationMdp
from LLC_Membranes.llclib import file_rw, topology, gromacs
from scipy.sparse import lil_matrix
import subprocess


def initialize():

    parser = argparse.ArgumentParser(description='Crosslink LLC structure')  # allow input from user

    parser.add_argument('-i', '--initial', default='initial.gro', help='Name of input file')
    parser.add_argument('-b', '--build_mon', default='NAcarb11Vd', type=str, help='Monomer the system is built with')
    parser.add_argument('-c', '--cutoff', default=5, type=float, help='Cutoff percentage for cross-linking. Bottom x %'
                                                              ' of the distribution of distances will be cross-linked ')
    parser.add_argument('-e', '--term_prob', default=5, type=float, help='Termination probability (%)')
    parser.add_argument('-d', '--dummy_name', default='dummies.gro', help='Name of initial .gro file after dummies are '
                                                                          'added')
    parser.add_argument('-itp', '--xlink_top_name', default='assembly.itp', help='Name of .itp topology file describing'
                        'assembly of cross-linked monomers.')
    parser.add_argument('-top', '--topname', default='topol.top', help='Name of topology file')
    parser.add_argument('-T', '--temperature', default=300, type=float, help='Temperature at which to run system (K)')
    parser.add_argument('--em_steps', default=5000, type=int, help='Maximum number of steps to run during energy'
                                                                  'minimization')
    parser.add_argument('-dt', '--dt', default=0.002, type=float, help='Time step for NVT simulations (ps)')
    parser.add_argument('-L', '--length', default=5, type=float, help='Length of NVT simulation (ps)')
    parser.add_argument('-ff', '--forcefield', default='gaff', help='Name of forcefield used during simulation')
    parser.add_argument('-mdp_em', '--mdp_em', default='em.mdp', help='Name of energy minimization .mdp file')
    parser.add_argument('-mdp_nvt', '--mdp_nvt', default='nvt.mdp', help='Name of energy minimization .mdp file')
    parser.add_argument('-res', '--residue', default='HII', help='Name of residue to be cross-linked')
    parser.add_argument('-resd', '--dummy_residue', default='HIId', help='Name of residue to be cross-linked with '
                                                                         'dummy atoms included in the topology')
    parser.add_argument('-density', '--density', default=95, type=float, help='Cross-link density (percent of x-links'
                        'that need to be cross-linked for procedure to terminate')
    parser.add_argument('-p', '--percent', default=1, type=float, help='Percent of eligible carbons that will bond')
    parser.add_argument('-rad_percent', default=20, type=float, help='Percent of radicals that react each iteration.'
                                                                     'Not stable above 50 % (?)')
    parser.add_argument('-out', '--output_gro', default='xlinked.gro', help='Name of final cross-linked structure')
    parser.add_argument('-rad_frac_term', default=20, type=float, help='Out of 100, how many radicals that will be '
                        'terminated on each iteration. Numbers greater than zero will decrease cross-link density')
    parser.add_argument('-stagnation', default=5, type=int, help='The number of iterations without generating a new'
                                                                 'cross-link before the algorithm forces termination')

    # parallelization options
    parser.add_argument('-mpi', '--parallelize', action="store_true", help='specify true if running with MPI or GPU')
    parser.add_argument('-np', '--nproc', default=4, type=int, help='Number of processess to run in parallel (number'
                                                                    'of GPUs on Bridges and Summit)')
    parser.add_argument('-dd', '--domain_decomposition', default=[2, 2, 1], help='xyz dimensions of domain '
                        'decomposition grid. This may need to be adjusted if there are issues with mdrun. The product'
                        'of these values should equal the number of processes (np)')

    # save options
    parser.add_argument('-s', '--save_intermediates', action="store_true", help='Save intermediate topology and energy'
                        'minimized structures in a folder called `intermediates`.')
    parser.add_argument('-sf', '--save_frequency', default=1, help='Number of iterations between saves')

    return parser


class Atoms(object):
    """ Give all atoms in the system properties including atom type, charge and mass

    :param res: Residue() object whose properties will be duplicated nres times
    :param nres: number of residues in full system

    :type res: Residue object
    :type nres: int
    """

    def __init__(self, res, nres):

        natoms = res.natoms
        self.name = np.zeros([nres*natoms], dtype=object)
        self.type = np.zeros([nres*natoms], dtype=object)
        self.charge = np.zeros([nres*natoms])
        self.mass = np.zeros([nres*natoms])

        for i in range(len(res.atom_info)):
            self.name[i::natoms] = str(res.atom_info[i][4])
            self.type[i::natoms] = str(res.atom_info[i][1])
            self.charge[i::natoms] = float(res.atom_info[i][6])
            self.mass[i::natoms] = float(res.atom_info[i][7].split(';')[0])


class Topology:

    def __init__(self):

        # create residue object for each residue in system
        self.residues = [topology.Residue(a, connectivity=True) for a in list(set([a.residue.name for a in
                                                                                   self.t.topology.atoms]))]

        # define which residue is going to be cross-linked
        self.res_ndx = [i.name for i in self.residues].index(self.xlink_residue_name)
        self.xlink_residue = self.residues[self.res_ndx]

        # number of atoms per residue
        self.nresidues = {}
        for r in self.residues:
            self.nresidues[r.name] = self.all_residues.count(r.name) // r.natoms

        # get all information about atoms
        self.xlink_residue_atoms = Atoms(self.xlink_residue, self.nresidues[self.xlink_residue_name])

        # get all bonds
        self.all_bonds = []
        natoms = self.xlink_residue.natoms
        bonds = self.xlink_residue.bonds
        for i in range(self.nresidues[self.xlink_residue_name]):
            for j in range(len(bonds)):
                self.all_bonds.append([int(bonds[j][0]) + i * natoms, int(bonds[j][1]) + i * natoms])

        # number all vsites
        nvsites = len(self.residues[self.res_ndx].virtual_sites)
        nres = self.nresidues[self.xlink_residue_name]
        self.vsites = np.zeros([nvsites*nres, 9], dtype=object)
        for i in range(nres):
            for k in range(nvsites):
                self.vsites[i * nvsites + k, :] = [float(j) for j in self.residues[self.res_ndx].virtual_sites[k]]
                self.vsites[i * nvsites + k, :5] += i * natoms

        # number all improper dihedrals
        nimpropers = len(self.xlink_residue.improper_dihedrals)
        self.impropers = np.zeros([nimpropers*nres, 5], dtype=int)
        for i in range(nres):
            for j in range(nimpropers):
                self.impropers[i * nimpropers + j, :] = [int(k) for k in self.xlink_residue.improper_dihedrals[j]]
                self.impropers[i * nimpropers + j, :4] += i * natoms

        # get pairs, angles, dihedrals
        self.pairs = None
        self.angles = None
        self.dihedrals = None
        self.adjacency_matrix = None

        self.define_topology()

        self.initiators = []  # indices of new initiator atoms which need to switched from virtual sites
        self.radicals = []  # indices of carbon atoms that are now radicals
        self.terminate = []

    def define_topology(self):

        self.create_adjacency_matrix()
        self.find_angles()
        self.find_dihedrals()
        self.find_pairs()

    def write_assembly_topology(self, out='assembly.itp', virtual_sites=True, vsite_atom_name='hc_d'):

        if not virtual_sites:
            count = 0
            new_numbers = {}  # everything will need to be renumbered

        with open(out, 'w') as f:

            f.write(';Topology for cross-linked assembly of %s residues\n' % self.xlink_residue_name)
            f.write('\n[ moleculetype ]\n')
            f.write(';name              nrexcl\n')
            f.write('%s                3\n' % self.xlink_residue_name)
            f.write('\n')
            f.write('[ atoms ]\n')
            f.write(';   nr  type  resi  res  atom  cgnr     charge      mass\n')
            atom_info = self.residues[self.res_ndx].atom_info
            natoms = self.residues[self.res_ndx].natoms
            nres = self.nresidues[self.xlink_residue_name]
            for i in range(nres):
                for j in range(len(atom_info)):
                    if virtual_sites:
                        f.write('{:>5d}{:>5s}{:>6s}{:>6s}{:>6s}{:>7s}{:>13f}{:>13f}'.format(int(atom_info[j][0]) + i*natoms,
                                self.xlink_residue_atoms.type[i*natoms + j], atom_info[j][2], atom_info[j][3],
                                atom_info[j][4], atom_info[j][5], self.xlink_residue_atoms.charge[i*natoms + j],
                                self.xlink_residue_atoms.mass[i * natoms + j]))
                        if (i * natoms + j) in self.radicals:
                            f.write('; *\n')
                        elif (i * natoms + j) in self.terminate:
                            f.write('; T\n')
                        else:
                            f.write('\n')
                    else:
                        if self.xlink_residue_atoms.type[i*natoms + j] != vsite_atom_name:
                            count += 1
                            f.write('{:>5d}{:>5s}{:>6s}{:>6s}{:>6s}{:>7s}{:>13f}{:>13f}\n'.format(count,
                                    self.xlink_residue_atoms.type[i*natoms + j], atom_info[j][2], atom_info[j][3],
                                    atom_info[j][4], atom_info[j][5], self.xlink_residue_atoms.charge[i*natoms + j],
                                    self.xlink_residue_atoms.mass[i * natoms + j]))
                            new_numbers[i*natoms + j + 1] = count

            if not virtual_sites:  # renumber atoms

                for i in range(len(self.all_bonds)):
                    self.all_bonds[i] = [new_numbers[x] for x in self.all_bonds[i]]
                for i in range(len(self.pairs)):
                    self.pairs[i] = [new_numbers[x] for x in self.pairs[i]]
                for i in range(len(self.angles)):
                    self.angles[i] = [new_numbers[x] for x in self.angles[i]]
                for i in range(len(self.dihedrals)):
                    self.dihedrals[i] = [new_numbers[x] for x in self.dihedrals[i]]
                for i in range(len(self.impropers)):
                    self.impropers[i, :-1] = [new_numbers[x] for x in self.impropers[i, :-1]]

            f.write('\n[ bonds ]\n')
            f.write(';   ai     aj funct\n')
            for i in range(len(self.all_bonds)):
                    f.write('{:>6d}{:>7d}  1\n'.format(self.all_bonds[i][0], self.all_bonds[i][1]))

            f.write('\n[ pairs ]\n')
            f.write(';   ai     aj    funct\n')
            for i in range(len(self.pairs)):
                f.write('{:>6d}{:>7d}{:>7d}\n'.format(self.pairs[i][0], self.pairs[i][1], 1))

            f.write('\n[ angles ]\n')
            f.write(';   ai      aj      ak    funct\n')
            for i in range(len(self.angles)):
                f.write('{:>6d}{:>7d}{:>7d}{:>7d}\n'.format(self.angles[i][0], self.angles[i][1], self.angles[i][2], 1))

            f.write('\n[ dihedrals ]\n')
            f.write(';   ai     aj    funct\n')
            for i in range(len(self.dihedrals)):
                f.write('{:>6d}{:>7d}{:>7d}{:>7d}{:>7d}\n'.format(self.dihedrals[i][0], self.dihedrals[i][1],
                                                                  self.dihedrals[i][2], self.dihedrals[i][3], 3))

            f.write('\n[ dihedrals ] ; impropers\n')
            f.write(';     i      j      k      l    func\n')
            #impropers = self.residues[self.res_ndx].improper_dihedrals

            for i in range(self.impropers.shape[0]):
                f.write('{:<6d}{:<7d}{:<7d}{:<7d}{:<7d}\n'.format(self.impropers[i, 0], self.impropers[i, 1],
                        self.impropers[i, 2], self.impropers[i, 3], self.impropers[i, 4]))

            if virtual_sites:
                f.write('\n[ virtual_sites4 ]\n')
                f.write(';Site  from                         funct   a          b          d\n')
                for i in range(self.vsites.shape[0]):
                    f.write('{:<8d}{:<6d}{:<6d}{:<6d}{:<8d}{:<8s}{:<11f}{:<11f}{:<11f}\n'.format(int(self.vsites[i, 0]),
                            int(self.vsites[i, 1]), int(self.vsites[i, 2]), int(self.vsites[i, 3]), int(self.vsites[i, 4]),
                            str(self.vsites[i, 5]), self.vsites[i, 6], self.vsites[i, 7], self.vsites[i, 8]))

    def create_adjacency_matrix(self):

        natoms = self.nresidues[self.xlink_residue_name] * self.residues[self.res_ndx].natoms
        #self.adjacency_matrix = np.zeros([natoms, natoms], dtype=bool)  # 27 GB --> 3 GB when changed dtype to bool
        self.adjacency_matrix = lil_matrix((natoms, natoms), dtype=int)  # 27 GB --> 56 bytes wow
        for i, bond in enumerate(self.all_bonds):
            self.adjacency_matrix[bond[0] - 1, bond[1] - 1] = 1
            self.adjacency_matrix[bond[1] - 1, bond[0] - 1] = 1

    def find_angles(self):

        self.angles = []
        for bond in self.all_bonds:
            # nonzero = np.nonzero(self.adjacency_matrix[bond[0] - 1])[0] + 1
            # nonzero = self.adjacency_matrix[bond[0] - 1].nonzero()[1] + 1
            nonzero = np.array(self.adjacency_matrix[bond[0] - 1].rows[0]) + 1
            for i in nonzero:
                if i < bond[1] and i not in bond:
                    self.angles.append([i, bond[0], bond[1]])

    def find_dihedrals(self):

        self.dihedrals = []
        for angle in self.angles:
            # nonzero = np.nonzero(self.adjacency_matrix[angle[0] - 1])[0] + 1
            #nonzero = self.adjacency_matrix[angle[0] - 1].nonzero()[1] + 1
            nonzero = np.array(self.adjacency_matrix[angle[0] - 1].rows[0]) + 1
            for i in nonzero:
                if i < angle[2] and i not in angle:
                    self.dihedrals.append([i, angle[0], angle[1], angle[2]])

    def find_pairs(self):

        self.pairs = []
        for i in self.dihedrals:
            self.pairs.append([i[0], i[3]])


class System(Topology):

    def __init__(self, initial_configuration, residue, dummy_residue, dummy_name='dummies.gro',
                 reaction_percentage=1, cutoff=0.6, radical_reaction_percentage=20, radical_termination_fraction=50):
        """ Initialize system for cross-linking and create a configuration with dummy atoms in appropriate locations

        :param initial_configuration: coordinate file to be cross-linked (.gro)
        :param residue: name of residue involved in cross-linking
        :param dummy_residue: file with same topology as residue in addition to dummy atoms
        :param dummy_name: name of coordinate file with dummy atoms
        :param reaction_percentage: percentage of atom pairs within the cutoff distance that will be chosen to bond
        :param cutoff: max distance two atoms can be away from each in order to still be considered from bonding
        :param radical_reaction_percentage: percentage of radicals that react on each iteration
        :param radical_termination_fraction: percentage of radicals that terminate each iteration

        :type initial_configuration: str
        :type residue: str
        :type dummy_residue: str
        :type dummy_name: str
        :type reaction_percentage: float
        :type cutoff: float
        :type radical_reaction_percentage: float
        :type radical_termination_fraction: float
        """

        add_dummies(initial_configuration, residue, dummy_residue,
                    out=dummy_name)  # add dummy atoms to the initial configuration

        self.reaction = XlinkReaction(residue)
        self.original_residue_name = residue
        self.xlink_residue_name = dummy_residue
        self.t = md.load(dummy_name)
        self.all_residues = [a.residue.name for a in self.t.topology.atoms]
        self.atom_names = [a.name for a in self.t.topology.atoms]

        super().__init__()  # inherit from Topology(). Topology() uses attributes defined above

        # self.n_xlink_residue = self.xlink_residue_atoms.name.size // self.xlink_residue.natoms

        self.chain1_atoms, self.chain2_atoms = [], []

        for i in self.reaction.scheme.monomer.bonds_with:
            self.chain1_atoms.append([a.index for a in self.t.topology.atoms if a.name in i[0]
                                      and a.residue.name == self.xlink_residue.name])
            self.chain2_atoms.append([a.index for a in self.t.topology.atoms if a.name in i[1]
                                      and a.residue.name == self.xlink_residue.name])

        self.bond_chain1, self.bond_chain2 = [], []
        self.reaction_type = []

        # Mostly informational variables
        self.nxlinks = 0
        # self.total_possible_xlinks = len(self.c1_atoms)
        # self.total_possible_terminated = len(self.c1_atoms) + len(self.c2_atoms)
        self.total_possible_xlinks = sum([len(i) for i in self.chain1_atoms])  # TODO: need to think about this more.
        self.total_possible_terminated = sum([len(i) for i in self.chain1_atoms]) + \
                                         sum([len(i) for i in self.chain2_atoms])
        self.iteration = 1
        self.terminated_radicals = []
        self.n_reacted_radicals = 0  # number of reacted radicals for a given iteration

        # Variables that control the reaction
        self.radical_reaction_percentage = radical_reaction_percentage / 100.  # percent of radicals that will react
        self.reaction_percentage = reaction_percentage / 100.
        self.cutoff = cutoff
        self.former_radicals = None
        self.rad_frac_term = radical_termination_fraction / 100.
        self.chain_numbers = []  # this is used to ensure that the same atoms don't react twice

        # Variables for topology modification
        self.rm_impropers = []

    def chain_number(self, ndx):
        """ Determine the chain number based on serial index """

        name = self.xlink_residue_atoms.name[ndx]
        mon = ndx // self.xlink_residue.natoms  # monomer number
        chain = self.reaction.scheme.determine_chains(name)[0]  # chain label

        return mon * self.reaction.scheme.monomer.nchains + self.reaction.scheme.monomer.chain_numbers[chain]

    def minimum_image_distances(self, list1, list2):
        """
        Find the minimum image distance between all pairs of atoms between two lists. They can be of different length.
        This is written for a cubic box, so transformations need to be applied if non-cubic box is used
        :param list1: numpy array or list
        :param list2: numpy array or list
        :return:
        """

        dim = self.t.unitcell_lengths[0]

        list1_len = len(list1)
        list2_len = len(list2)

        diff = np.zeros([list1_len * list2_len, 3])
        for i in range(list1_len):
            diff[i * list2_len: (i + 1) * list2_len] = self.t.xyz[0, list1[i], :] - self.t.xyz[0, list2, :]

        diff = np.where(diff > 0.5 * dim, diff - dim, diff)
        diff = np.where(diff < -0.5 * dim, diff + dim, diff)

        return diff

    def generate_ordered_distances(self, list1, list2):
        """
        Generate a list ordered sequentially based on pairwise distances between atom indices listed in list1 and list2
        :param list1: list or numpy array
        :param list2: list or numpy array
        :param cut: maximum allowable pairwise distance to still be considered eligible for bonding (nm)
        :return: 2 same-length numpy arrays ordered so each index is paired with the same index with the other array.
        Pairwise distances are in increasing order.
        """

        diff = self.minimum_image_distances(list1, list2)
        d = np.linalg.norm(diff, axis=1)

        d_ndx = np.argsort(d)

        # These lists need post-modification to account for adjacent atoms
        eligible_list1 = (d_ndx // len(list2))
        eligible_list2 = (d_ndx % len(list2))

        eligible_list1 = np.array(list1)[eligible_list1]
        eligible_list2 = np.array(list2)[eligible_list2]

        cut_ndx = np.argmin(np.abs(d[d_ndx] - self.cutoff))  # index where cut is exceeded

        return eligible_list1[:cut_ndx], eligible_list2[:cut_ndx]

    def select_eligible_carbons(self):
        """
        - Calculate the distance between atoms of interest.
        - Find out which pairs meet the distance cutoff criteria and make realistic bonds
        - Determine type of reaction
        """

        # calculate minimum image distance between c1 and c2 atoms
        if not np.isclose(self.t.unitcell_angles[:, 2].mean(), 90., rtol=1e-01):
            # If unit cell is non-cubic, transform it to cubic -- this is not generalized. Works for HII unit cell
            print('hellurr')
            # shift so box is temporarily cubic. The min image distances are just an approximation using this approach
            theta = self.t.unitcell_angles[:, 2].mean() * (np.pi / 180.)

            self.t.xyz[..., 1] /= np.sin(theta)
            self.t.xyz[..., 0] = self.t.xyz[..., 0] - self.t.xyz[..., 1]*np.cos(theta)

        if self.radicals:
            eligible_c_rad, eligible_rad = [], []
            for i, c in enumerate(self.chain1_atoms):
                d = self.generate_ordered_distances(c + self.chain2_atoms[i], self.radicals)
                eligible_c_rad += d[0].tolist()  # list of numpy arrays
                eligible_rad += d[1].tolist()

            bonds = self.bond_filter(np.array(eligible_c_rad), np.array(eligible_rad), radical=True)
            self.bond_chain1 += bonds[0]
            self.bond_chain2 += bonds[1]
            for i, x in enumerate(self.bond_chain1):
                print(x, self.bond_chain2[i])
            print('----------------------------------')

            # eligible_c1_rad, eligible_rad = self.generate_ordered_distances(self.chain1_atoms, self.radicals)
            # bonds = self.bond_filter(eligible_c1_rad, eligible_rad, radical=True)
            # self.bond_chain1 += bonds[0]
            # self.bond_chain2 += bonds[1]
            # for i in range(len(bonds[0])):
            #     self.reaction_type.append('radical_c2')

        eligible_c1, eligible_c2 = [], []
        for i, c in enumerate(self.chain1_atoms):

            d = self.generate_ordered_distances(c, self.chain2_atoms[i])
            eligible_c1 += d[0].tolist()  # list of numpy arrays
            eligible_c2 += d[1].tolist()

        bonds = self.bond_filter(np.array(eligible_c1), np.array(eligible_c2))
        self.bond_chain1 += bonds[0]
        self.bond_chain2 += bonds[1]

        for i, x in enumerate(self.bond_chain1):
            print(x, self.bond_chain2[i], self.reaction_type[i])

        # for i in range(len(bonds[0])):
        #     self.reaction_type.append('head2tail')

        # eligible_c1, eligible_c2 = self.generate_ordered_distances(self.c1_atoms, self.c2_atoms)
        # nreact = int(self.reaction_percentage * len(eligible_c1))
        # bonds = self.bond_filter(eligible_c1, eligible_c2, nreact)
        # self.bond_c1 += bonds[0]
        # self.bond_c2 += bonds[1]
        # for i in range(len(bonds[0])):
        #     self.reaction_type.append('head2tail')

    def bond_filter(self, eligible_chain1, eligible_chain2, radical=False):
        """ Choose which carbons to bond based on proximity. In the process filter out bonds that should not form --
        two bonds to the same carbon, bonds between carbons on same chain, simultaneous bonds on same chain.

        :param eligible_chain1: chain1 atoms to bond. Formatted as a numpy array for efficient slicing
        :param eligible_chain2: chain2 atoms to bond. This list determines the number of bonds based on the percent
        :param radical: True if this is a radical reaction. A different percentage of eligible chains will react

        :type eligible_chain1: np.ndarray
        :type eligible_chain2: np.ndarray
        :type radical: bool
        """

        # Determine reaction types and eliminate potential reactions that don't actually happen. So far this only
        # happens with radicals
        rxn_types = np.zeros([len(eligible_chain1)], dtype=object)
        keep = []
        for i, x in enumerate(eligible_chain1):
            c1_name = self.atom_names[x]
            c2_name = self.atom_names[eligible_chain2[i]]
            rxn_types[i] = self.reaction.scheme.determine_reaction_type(c1_name, c2_name, radical=radical)
            if rxn_types[i]:
                keep.append(i)

        eligible_chain1 = eligible_chain1[keep]
        eligible_chain2 = eligible_chain2[keep]
        rxn_types = rxn_types[keep]

        # prevent same-chain carbon atoms from bonding
        keep = []
        for i, x in enumerate(eligible_chain1):
            if not self.check_pair(x, eligible_chain2[i]):
                keep.append(i)

        eligible_chain1 = eligible_chain1[keep]
        eligible_chain2 = eligible_chain2[keep]
        rxn_types = rxn_types[keep]

        # narrow list to bottom percent of distances
        # redetermine nreact based on number of eligible carbons _after_ above loop
        percent = self.reaction_percentage
        if radical:
            percent = self.radical_reaction_percentage
        nreact = int(percent * len(eligible_chain2))  # total number of reactions

        reactions = self.reaction.scheme.reaction_weights
        if radical:
            reactions = self.reaction.scheme.radical_reaction_weights

        n_rxn_per_type = {k: int(reactions[k]*nreact) for k in list(reactions.keys())}  # number reactions of each type
        n_rxn_count = {k: 0 for k in list(reactions.keys())}  # counter for number of reactions of each type so far

        bond_chain1 = []
        bond_chain2 = []
        for i, r in enumerate(rxn_types):
            if n_rxn_count[r] < n_rxn_per_type[r]:
                chain1 = self.chain_number(eligible_chain1[i])
                chain2 = self.chain_number(eligible_chain2[i])
                if chain1 == chain2:  # this shouldn't happen, but this will let me know if it does.
                    print('chain1 equals chain2!')
                if chain1 not in self.chain_numbers and chain2 not in self.chain_numbers:  # make sure its always new independent chains reacting
                    bond_chain1.append(eligible_chain1[i] + 1)  # convert to serial
                    bond_chain2.append(eligible_chain2[i] + 1)
                    n_rxn_count[r] += 1
                    self.reaction_type.append(r)
                    self.chain_numbers += [chain1, chain2]

        return bond_chain1, bond_chain2

    def check_pair(self, i, j):
        """ Determine if two potential bonding carbon atoms are on the same chain

        :return: True if pair is on same chain, False of pair is part of different chains
        """

        # check if carbons are a part of the same monomer
        if i // self.xlink_residue.natoms == j // self.xlink_residue.natoms:
            # check if carbons are part of the same chain
            chain1, chain2 = self.reaction.scheme.determine_chains([self.atom_names[i], self.atom_names[j]])
            if chain1 == chain2:
                return True
            else:
                return False
        else:
            return False

    def bond(self, quench=True):
        """ Create all new cross-links and terminate some of the radicals that are formed

        :param quench: if True, do not try to bond anything and just terminate all remaining radicals

        :type quench: bool
        """

        if not quench:

            for i, x in enumerate(self.bond_chain1):

                c1_name = self.xlink_residue_atoms.name[x - 1]  # might be better to use self.atom_names
                c2_name = self.xlink_residue_atoms.name[self.bond_chain2[i] - 1]

                # # add chain1--chain2 bond
                # self.all_bonds.append([x, self.bond_chain2[i]])

                self.react(self.reaction_type[i], {c1_name: x, c2_name: self.bond_chain2[i]})

        self.terminate_radicals(quench=quench)  # this needs to go here since termination requires an improper to be removed

        rm_improper_ndx = []
        for i, p in enumerate(self.impropers):
            if p[:4].tolist() in self.rm_impropers:
                rm_improper_ndx.append(i)

        self.impropers = np.delete(self.impropers, rm_improper_ndx, axis=0)

        self.nxlinks += len(self.bond_chain1)  # update total number of cross-links in system

        # self.identify_terminated(rad_term_frac=self.rad_frac_term)  # this will become it's own reaction scheme
        self.remove_virtual_sites()  # could go in self.react and get rid of self.initiators

        self.define_topology()  # make new angles, pairs and dihedrals lists

    def react(self, reaction_type, atoms):
        """ React according the reaction type and scheme

        :param reaction_type: type of reaction
        :param atoms: dictionary of atom names and their corresponding serial indices

        :type reaction_type: str
        :type atoms: dict
        """

        types, bonds, radicals, rm_impropers, terminate = self.reaction.scheme.react(reaction_type, atoms)
        # print(types)
        # print(bonds)
        # print(radicals)
        # print(rm_impropers)
        # print(terminate)

        # add bonds between dummy atoms (made real) and carbon atoms
        for b in bonds:
            self.all_bonds.append(b[:2])
            if b[2] == 'dummy':
                self.xlink_residue_atoms.mass[b[1] - 1] = self.reaction.scheme.monomer.dummy_mass  # dummy always 2nd entry in b
                self.initiators.append(b[1])  # The hydrogen dummy atoms are our 'initiators'.

        # change atom types of newly bonded carbon and former dummy atoms
        for k in types.keys():
            for j in types[k].keys():
                self.xlink_residue_atoms.type[j - 1] = types[k][j]

        # add new radical to list of radicals
        for r in radicals:
            self.radicals.append(r - 1)

        # remove improper dihedrals where hybridization changed
        for p in rm_impropers:
            self.rm_impropers.append(p)

        # terminate sp3 things
        for t in terminate:
            self.terminate.append(t - 1)

        # remove radical if this is a radical reaction -- this will need to be improved to be more general. 
        # It could be a return in the reaction schemes
        if 'radical' in reaction_type:
            ndx = list(atoms.values())[1] - 1
            self.radicals.remove(ndx)
            self.terminated_radicals.append(ndx)

        if 'terminate' in reaction_type:
            ndx = terminate[0] - 1
            self.radicals.remove(terminate[0] - 1)
            self.terminated_radicals.append(ndx)

    def terminate_radicals(self, quench=False):
        """ Choose which remaining radicals to terminate and then terminate them. All reacted radicals should have been
        removed before running this.

        :param quench: if True, terminate all remaininng radicals

        :type quench: bool
        """

        if quench:
            terminate = self.radicals
        else:
            nterm = int(self.rad_frac_term * len(self.radicals))
            terminate = np.random.choice(self.radicals, size=nterm, replace=False)

        for t in terminate:

            name = self.xlink_residue_atoms.name[t]  # self.radicals indexing starts at 0

            self.react('terminate', {name: t + 1})

    def remove_virtual_sites(self):

        # remove virtual sites for atoms that became real
        rm_vsite = []
        for i in range(self.vsites.shape[0]):
            if int(self.vsites[i, 0]) in self.initiators:
                rm_vsite.append(i)

        self.vsites = np.delete(self.vsites, rm_vsite, axis=0)  # delete vsite entries

    def reload_coordinates(self, name):

        self.t = md.load(name)

    def update_lists(self):

        # reacted chain1 and chain2 atoms are no longer eligible for bonding
        # self.chain1_atoms is numbered by index, while self.bond_chain1 is number serially
        for i, c in enumerate(self.chain1_atoms):
            self.chain1_atoms[i] = [x for x in c if x not in self.terminate and x not in self.radicals]
            self.chain2_atoms[i] = [x for x in self.chain2_atoms[i] if x not in self.terminate and x not in
                                    self.radicals]

        # re-initialize bonding carbon lists
        self.bond_chain1 = []
        self.bond_chain2 = []
        self.reaction_type = []
        self.chain_numbers = []

        self.initiators = []  # bad name for this. It's really just dummy atoms that become real
        self.rm_impropers = []
        self.terminated_radicals = []

    def update_log(self):

        write_mode = 'a'
        if self.iteration == 1:
            write_mode = 'w'

        with open('xlink.log', '%s' % write_mode) as f:

            f.write('Iteration Number: %d \n' % self.iteration)
            f.write('Number of new cross-links: %d\n' % len(self.bond_chain1))
            f.write('Total system cross-links: %d\n' % self.nxlinks)
            f.write('Cross-link density: %.2f %%\n' % (100. * (self.nxlinks / self.total_possible_xlinks)))
            f.write('+-----------------------------------+\n')
            f.write('| Serial indices of new cross-links |\n')
            f.write('+------+------+---------------------+\n')
            f.write('|{:^6}|{:^6}|{:^21}|\n'.format('chain1', 'chain2', 'Reaction Type'))
            f.write('+------+------+\n')
            for i, x in enumerate(self.bond_chain1):
                f.write('{:>5d} -- {:<5d}  {}\n'.format(x, self.bond_chain2[i], self.reaction_type[i]))
            f.write('Total radicals terminated: %.2d\n' % len(self.terminated_radicals))
            f.write('+---------------------------------------+\n')
            f.write('| Serial indices of terminated radicals |\n')
            f.write('+------+------+-------------------------+\n')
            for i, x in enumerate(self.terminated_radicals):
                if i != 0 and i % 10 == 0 or i + 1 == len(self.terminated_radicals):
                    f.write('{:>6d}\n'.format(x))
                else:
                    f.write('{:>6d}'.format(x))
            f.write('Total radicals left in system: %d\n' % len(self.radicals))
            f.write('Percent terminated: %.2f %%\n' % (100 * len(self.terminate) / (self.total_possible_terminated)))
            f.write('%s\n' % (80*'-'))

    def cleanup(self, out='xlinked.gro', dummy_atom_name='hc_d'):
        """ Remove virtual sites in .gro and topology.

        :param out: name of output .gro file to be written
        :param dummy_atom_name: atom type of dummy atom
        :return: A .gro file and an .itp topology file without dummy atoms
        """

        keep = [i for i, x in enumerate(self.xlink_residue_atoms.type) if x != dummy_atom_name]
        keep += [a.index for a in self.t.topology.atoms if a.residue.name != self.original_residue_name]

        ids = np.array([a.name for a in self.t.topology.atoms], dtype=object)[keep]
        res = np.array([a.residue.name for a in self.t.topology.atoms], dtype=object)[keep]

        # mdtraj workaround because SOL gets renamed to HOH, OW1 to O, HW1 to H1 and HW2 to H2
        for i, x in enumerate(ids):
            if res[i] == 'HOH':
                if x == 'O':
                    ids[i] = 'OW'
                if x == 'H1':
                    ids[i] = 'HW1'
                if x == 'H2':
                    ids[i] = 'HW2'
                res[i] = 'SOL'

        ucell = self.t.unitcell_vectors[-1, ...]

        file_rw.write_gro_pos(self.t.xyz[-1, keep, :], out, ids=ids, res=res, ucell=ucell)

        self.write_assembly_topology(virtual_sites=False, vsite_atom_name=dummy_atom_name)


def save_intermediates(files, iteration):
    """ For a given iteration, save intermediate files in a folder called `intermediates`

    :param files: list of file names to save
    :param iteration: iteration number used for naming

    :type file: list of str
    :type iteration: int
    """

    if not os.path.isdir("./intermediates"):
        os.mkdir('intermediates')

    for f in files:
        cp = "cp %s intermediates/%s_%d" % (f, f, iteration)
        p = subprocess.Popen(cp.split())
        p.wait()


def crosslink(params):
    """ The main cross-linking algorithm

    :param params: a dictionary of parameters controlling the cross-linking reaction. There is a dictionary key for
    every argument in xlink.initialize()

    :type params: dict
    """

    start = time.time()
    initial_message = '# Cross-linking %s #' % params['initial']
    print('#' * len(initial_message))
    print(initial_message)
    print('#' * len(initial_message))

    # if needed
    dt_short = 0.0005  # simulation time step (ps) (unstable if timestep is too large)
    length_short = .025  # (integer for now, so this is as low as it goes) length of short nvt simulation in picoseconds

    # get the system ready for cross-linking
    # Add dummy atoms create topology for entire assembly
    print('Initializing system and adding dummy atoms...', end='', flush=True)
    sys = System(params['initial'], params['residue'], params['dummy_residue'], dummy_name=params['dummy_name'],
                 radical_reaction_percentage=params['rad_percent'], reaction_percentage=params['percent'],
                 radical_termination_fraction=params['rad_frac_term'])
    print('Done!\n%s created' % params['dummy_name'])

    print('Generating input files: %s, %s, %s and %s...' % (params['xlink_top_name'], params['topname'],
                                                            params['mdp_em'], params['mdp_nvt']), end='', flush=True)
    sys.write_assembly_topology()

    # generate input files that will be used throughtout cross-linking process
    full_top = SystemTopology(params['dummy_name'], ff=params['forcefield'], xlink=True,
                              xlinked_top_name=params['xlink_top_name'])  # topology with other residues and forcefield
    full_top.write_top(name=params['topname'])
    mdpshort = SimulationMdp(params['dummy_name'], T=params['temperature'], em_steps=params['em_steps'],
                             time_step=params['dt'], length=params['length'], p_coupling='semiisotropic',
                             barostat='berendsen', genvel='yes', restraints=False, xlink=True, bcc=False)
    mdpshort.write_em_mdp(out='%s' % params['mdp_em'].split('.')[0])
    mdpshort.write_nvt_mdp(out='%s' % params['mdp_nvt'].split('.')[0])

    # TODO: I'm not sure the short time step NVT simulation is used (or needed)
    mdpshort = SimulationMdp(params['dummy_name'], T=params['temperature'], em_steps=params['em_steps'],
                             time_step=dt_short, length=length_short, p_coupling='semiisotropic', barostat='berendsen',
                             genvel='yes', restraints=False, xlink=True, bcc=False)
    mdpshort.write_nvt_mdp(out='nvt_short')
    print('Done!')

    # energy minimize starting configuration to get dummies in the right place
    print('Energy minimizing %s...' % params['dummy_name'], end='', flush=True)
    gromacs.simulate(params['mdp_em'], params['topname'], params['dummy_name'],
                     'em_%s' % params['dummy_name'].split('.')[0], mpi=params['parallelize'],
                     nprocesses=params['nproc'], dd=params['domain_decomposition'])
    print('Done!')

    # Rest of iterations
    stagnated_iterations = 0  # number of iterations without forming a cross-link
    while (len(sys.terminate) / sys.total_possible_terminated) < (params['density'] / 100.) and stagnated_iterations < \
            params['stagnation']:

        #### FOR TESTING ####
        # sys.iteration = 2
        # lists = np.load('radicals.npz')
        # sys.radicals = list(lists['radicals'])
        # sys.terminate = list(lists['term'])
        #####################

        print('-' * 80)
        print('Iteration: %d' % sys.iteration)

        if sys.iteration > 1:
            # update coordinates
            sys.reload_coordinates('nvt.gro')
            sys.update_lists()

        print('Choosing atoms and cross-linking them...', end='', flush=True)
        sys.select_eligible_carbons()
        sys.bond()
        sys.write_assembly_topology()  # re-write assembly topology
        print('Done!')
        print('Total new bonds: %d' % len(sys.bond_chain1))
        np.savez_compressed('radicals', radicals=sys.radicals, term=sys.terminate)

        # energy minimize then run short NVT simulation
        print('Energy minimizing new cross-links...', end='', flush=True)
        if sys.iteration == 1:
            gromacs.simulate(params['mdp_em'], params['topname'], 'em_%s' % params['dummy_name'], 'em',
                             mpi=params['parallelize'], nprocesses=params['nproc'], dd=params['domain_decomposition'])
        else:
            gromacs.simulate(params['mdp_em'], params['topname'], 'nvt.gro', 'em', mpi=params['parallelize'],
                             nprocesses=params['nproc'], dd=params['domain_decomposition'])
        print('Done!')

        print('Running %.1f ps NVT simulation...' % params['length'], end='', flush=True)
        gromacs.simulate(params['mdp_nvt'], params['topname'], 'em.gro', 'nvt', mpi=params['parallelize'],
                         nprocesses=params['nproc'], dd=params['domain_decomposition'])
        print('Done!')

        print('\nTotal new cross-links: %d' % len(sys.bond_chain1))
        print('Total system cross-links: %d' % sys.nxlinks)
        print('Cross-link density: %.2f %%' % (100. * (sys.nxlinks / sys.total_possible_xlinks)))
        print('Total radicals terminated: %.2d' % len(sys.terminated_radicals))
        print('Total radicals left in system: %d' % len(sys.radicals))
        print('Percent terminated: %.2f' % (100 * len(sys.terminate) / sys.total_possible_terminated))

        if params['save_intermediates'] and sys.iteration % params['save_frequency'] == 0:
            save_intermediates(['em.gro', 'nvt.gro', 'assembly.itp'], sys.iteration)

        if len(sys.bond_chain1) > 0:
            stagnated_iterations = 0
        else:
            stagnated_iterations += 1

        # update log
        sys.update_log()
        sys.iteration += 1

    print('Quenching reaction...', end='', flush=True)
    sys.reload_coordinates('nvt.gro')
    sys.update_lists()
    sys.bond(quench=True)  # terminate all radicals
    print('Done!')

    print('Removing remaining dummy atoms from .gro and .itp files...', end='', flush=True)
    sys.cleanup(out=params['output_gro'])
    print('Done!')

    print('Energy minimizing %s...' % params['output_gro'], end='', flush=True)
    gromacs.simulate(params['mdp_em'], params['topname'], params['output_gro'], params['output_gro'].split('.')[0],
                     mpi=params['parallelize'], nprocesses=params['nproc'])
    print('Done!')

    print('Finished cross-linking in %.2f seconds' % (time.time() - start))


if __name__ == "__main__":

    args = initialize().parse_args()

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    crosslink(vars(args))
