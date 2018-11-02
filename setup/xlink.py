#!/usr/bin/env python

import mdtraj as md
import argparse
import time
import os
import numpy as np
import subprocess
from LLC_Membranes.setup.add_dummies import add_dummies
from LLC_Membranes.setup.gentop import SystemTopology
from LLC_Membranes.setup.genmdp import SimulationMdp
from LLC_Membranes.llclib import file_rw
from scipy.sparse import lil_matrix


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
                                                                     'Not stable above 50 %')
    parser.add_argument('-out', '--output_gro', default='xlinked.gro', help='Name of final cross-linked structure')
    parser.add_argument('-rad_frac_term', default=50, type=float, help='Out of 100, how many radicals that will be terminated'
                                                                        'on each iteration.')
    parser.add_argument('-stagnation', default=5, type=int, help='The number of iterations without generating a new'
                                                                 'cross-link before the algorithm forces termination')

    # parallelization options
    parser.add_argument('-mpi', '--parallelize', action="store_true", help='specify true if running with MPI or GPU')
    parser.add_argument('-np', '--nproc', default=4, type=int, help='Number of processess to run in parallel (number'
                                                                    'of GPUs on Bridges and Summit)')

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
        self.type = np.zeros([nres*natoms], dtype=object)
        self.charge = np.zeros([nres*natoms])
        self.mass = np.zeros([nres*natoms])

        for i in range(len(res.atom_info)):
            self.type[i::natoms] = str(res.atom_info[i][1])
            self.charge[i::natoms] = float(res.atom_info[i][6])
            self.mass[i::natoms] = float(res.atom_info[i][7])


class Residue(object):
    """ Get topology information about a given residue

    :param: name: name of residue. name.itp and name.gro should be stored in llcsim/top/topologies

    :type: name: str
    """

    def __init__(self, name):

        self.name = name
        self.is_ion = False
        self.check_ion()
        self.virtual_sites = None
        self.c1_atoms = None
        self.c2_atoms = None
        self.c1_atom_indices = None
        self.c2_atom_indices = None
        self.carbonyl_oxygen_indices = None

        if not self.is_ion:

            self.itp = []
            with open('%s/../top/topologies/%s.itp' % (location, name), 'r') as f:
                for line in f:
                    self.itp.append(line)

            self.bonds = self.get_bonds()

            self.atom_info = self.get_atoms()
            self.natoms = len(self.atom_info)

            self.virtual_sites = self.get_vsites()
            self.improper_dihedrals = self.get_improper_dihedrals()

        else:
            self.natoms = 1

    def check_ion(self):
        """ Check whether the given residue is an ion. Return True its name is stored in ions.txt located in
        llcsim/top/topologies

        """

        with open('%s/../top/topologies/ions.txt' % location) as f:
            ions = []
            for line in f:
                if line[0] != '#':
                    ions.append(str.strip(line))

        if self.name in ions:
            self.is_ion = True

    def get_atoms(self):
        """ Store all information inside the "[ atoms ]" section of name.itp

        """

        atoms_index = 0
        while self.itp[atoms_index].count('[ atoms ]') == 0:
            atoms_index += 1
        atoms_index += 1

        while self.itp[atoms_index].split()[0] == ';':
            atoms_index += 1

        atoms = []
        while self.itp[atoms_index] != '\n':
            atoms.append(self.itp[atoms_index].split())
            atoms_index += 1

        return atoms

    def get_bonds(self):
        """ Store all information in the "[ bonds ]" section of name.itp

        """

        bonds_index = 0
        while self.itp[bonds_index].count('[ bonds ]') == 0:
            bonds_index += 1
        bonds_index += 1

        while self.itp[bonds_index].split()[0] == ';':
            bonds_index += 1

        bonds = []
        while self.itp[bonds_index] != '\n':
            bonds.append([int(self.itp[bonds_index].split()[i]) for i in range(2)])
            bonds_index += 1

        return bonds

    def get_improper_dihedrals(self):
        """ Store all information in the "[ dihedrals ] ; impropers" section of name.itp

        """

        imp_ndx = 0

        while self.itp[imp_ndx].count('[ dihedrals ] ; impropers') == 0:
            imp_ndx += 1
            if imp_ndx >= len(self.itp):
                break

        impropers = None
        if imp_ndx < len(self.itp):
            impropers = []
            imp_ndx += 1
            while self.itp[imp_ndx][0] == ';':
                imp_ndx += 1
            while imp_ndx < len(self.itp) and self.itp[imp_ndx] != '\n':
                impropers.append(self.itp[imp_ndx].split())
                imp_ndx += 1

        return impropers

    def get_vsites(self):
        """ Store all information in the "[ virtual_sites ]" section of name.itp

        """

        vsite_index = 0
        while self.itp[vsite_index].count('[ virtual_sites4 ]') == 0:
            vsite_index += 1
            if vsite_index >= len(self.itp):
                break

        vsites = None
        if vsite_index < len(self.itp):
            vsites = []
            vsite_index += 1
            while self.itp[vsite_index][0] == ';':
                vsite_index += 1
            while vsite_index < len(self.itp):
                vsites.append(self.itp[vsite_index].split())
                vsite_index += 1

        return vsites

    def get_reactive_carbons(self):
        """ Read annotations in name.gro to determine which

        """

        self.c1_atoms = []
        self.c2_atoms = []
        self.c1_atom_indices = []
        self.c2_atom_indices = []

        with open('%s/../top/topologies/%s.gro' % (location, self.name), 'r') as f:

            for line in f:
                if line.count(';') != 0:
                    if line.split(';')[1].count('C1') > 0:
                        self.c1_atoms.append(line[10:15].strip())
                        self.c1_atom_indices.append(int(line[15:20]))
                    if line.split(';')[1].count('C2') > 0:
                        self.c2_atoms.append(line[10:15].strip())
                        self.c2_atom_indices.append(int(line[15:20]))

    def get_carbonyl_oxygens(self):

        self.carbonyl_oxygen_indices = []

        with open('%s/../top/topologies/%s.gro' % (location, self.name), 'r') as f:

            for line in f:
                if line.count(';') != 0:
                    if line.split(';')[1].count('O') > 0:
                        self.carbonyl_oxygen_indices.append(int(line[15:20]))

        return self.carbonyl_oxygen_indices


class Topology():

    def __init__(self):

        # create residue object for each residue in system
        self.residues = [Residue(a) for a in list(set([a.residue.name for a in self.t.topology.atoms]))]

        # define which residue is going to be cross-linked
        self.res_ndx = [i.name for i in self.residues].index(self.xlink_residue_name)
        self.xlink_residue = self.residues[self.res_ndx]
        self.xlink_residue.get_reactive_carbons()  # get names of carbon atoms that will bond

        # number of atoms per residue
        self.nresidues = {}
        for r in self.residues:
            self.nresidues[r.name] = self.all_residues.count(r.name) // r.natoms

        # get all information about atoms
        self.xlink_residue_atoms = Atoms(self.xlink_residue, self.nresidues[self.xlink_residue_name])

        # get all bonds
        self.all_bonds = []
        natoms = self.residues[self.res_ndx].natoms
        bonds = self.residues[self.res_ndx].bonds
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
                        elif (i * natoms + j + 1) in self.terminate:
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

        add_dummies(md.load(initial_configuration), residue, dummy_residue,
                    out=dummy_name)  # add dummy atoms to the initial configuration

        self.original_residue_name = residue
        self.xlink_residue_name = dummy_residue
        self.t = md.load(dummy_name)
        self.all_residues = [a.residue.name for a in self.t.topology.atoms]

        super().__init__()  # inherit from Topology()

        self.c1_atoms = [a.index for a in self.t.topology.atoms if a.name in self.xlink_residue.c1_atoms]
        self.c2_atoms = [a.index for a in self.t.topology.atoms if a.name in self.xlink_residue.c2_atoms]

        self.bond_c1 = []
        self.bond_c2 = []

        # Mostly informational variables
        self.nxlinks = 0
        self.total_possible_xlinks = len(self.c1_atoms)
        self.total_possible_terminated = len(self.c1_atoms) + len(self.c2_atoms)
        self.iteration = 1
        self.terminated_radicals = []
        self.n_reacted_radicals = 0  # number of reacted radicals for a given iteration

        # Variables that control the reaction
        self.radical_reaction_percentage = radical_reaction_percentage / 100.  # percent of radicals that will react
        self.reaction_percentage = reaction_percentage / 100.
        self.cutoff = cutoff
        self.former_radicals = None
        self.rad_frac_term = radical_termination_fraction / 100.

    def simulate(self, configuration, mdp='em.mdp', top='topol.top', out='em', parallel=True, np=4):
        """
        Energy minimize a configuration using existing .mdp files
        """

        p1 = subprocess.Popen(
            ["gmx", "grompp", "-p", "%s" % top, "-f", "%s" % mdp, "-o", "%s" % out, "-c", "%s" % configuration],
            stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # generate atomic level input file
        p1.wait()

        if parallel:
            p2 = subprocess.Popen(["mpirun", "-np", "%s" % np, "gmx_mpi", "mdrun", "-deffnm", "%s" % out],
                                  stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        else:
            p2 = subprocess.Popen(["gmx", "mdrun", "-deffnm", "%s" % out], stdout=open(os.devnull, 'w'),
                                  stderr=subprocess.STDOUT)  # run energy minimization
        p2.wait()

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

        # These lists may need post-modification to account for adjacent atoms
        eligible_list1 = (d_ndx // len(list2))
        eligible_list2 = (d_ndx % len(list2))

        eligible_list1 = np.array(list1)[eligible_list1]
        eligible_list2 = np.array(list2)[eligible_list2]

        cut_ndx = np.argmin(np.abs(d[d_ndx] - self.cutoff))  # index where cut is exceeded

        return eligible_list1[:cut_ndx], eligible_list2[:cut_ndx]

    def select_eligible_carbons(self):
        """
        This whole function could be written better since generating radical list is same code as bonded carbon list
        - Calculate the distance between atoms of interest. Find out which pairs meet the distance cutoff criteria
        - Head-to-tail addition dominates due to sterics and radical stabilization
        :param percent: the bottom percentage of the distribution of distances (below cut) that will be chosen for cross-linking
        :param cut: any distance beyond 'cut' will not be accounted for when calculating percent
        :param percent_rad: percentage of eligible radicals that will react
        """

        # calculate minimum image distance between c1 and c2 atoms

        # shift so box is cubic temporarily. The minimum image distances are just an approximation using this approach
        theta = np.pi / 3

        self.t.xyz[..., 1] /= np.sin(theta)
        self.t.xyz[..., 0] = self.t.xyz[..., 0] - self.t.xyz[..., 1]*np.cos(theta)

        if self.radicals:

            eligible_c1_rad, eligible_rad = self.generate_ordered_distances(self.c1_atoms, self.radicals)

            nreact = int(self.radical_reaction_percentage * len(self.radicals))

            for i, x in enumerate(eligible_c1_rad):
                if len(self.bond_c1) < nreact and (x + 1) not in self.bond_c1 and (eligible_rad[i] + 1) not in self.bond_c2:
                    self.bond_c1.append(x + 1)
                    self.bond_c2.append(eligible_rad[i] + 1)
                if len(self.bond_c1) >= nreact:
                    break

        eligible_c1, eligible_c2 = self.generate_ordered_distances(self.c1_atoms, self.c2_atoms)

        # prevent adjacent carbon atoms from bonding
        exclude = 0
        while eligible_c1[exclude] == (eligible_c2[exclude] + 1):  # This 1 is hard-coded. Could be made more general
            exclude += 1

        eligible_c1 = eligible_c1[exclude:]
        eligible_c2 = eligible_c2[exclude:]

        # narrow list to bottom percent of distances
        stop = int(self.reaction_percentage * len(eligible_c1))
        self.bond_c1 += (eligible_c1[:stop] + 1).tolist()  # convert indices to serial
        self.bond_c2 += (eligible_c2[:stop] + 1).tolist()

        # do not allow a carbon to bond twice (prioritizes radical reactions)
        keep_c1 = np.sort(np.unique(self.bond_c1, return_index=True)[1])
        self.bond_c1 = np.array(self.bond_c1)[keep_c1]
        self.bond_c2 = np.array(self.bond_c2)[keep_c1]

        keep_c2 = np.sort(np.unique(self.bond_c2, return_index=True)[1])
        self.bond_c1 = self.bond_c1[keep_c2]
        self.bond_c2 = self.bond_c2[keep_c2]

        # ensure that a tail is not donating and recieving bonds at the same time
        delete = []
        for i, x in enumerate(self.bond_c1[::-1]):  # go backwards to maintain radical prioritization
            if x - 1 in self.bond_c2:
                delete.append(len(self.bond_c1) - i - 1)

        self.bond_c1 = np.delete(self.bond_c1, delete)
        self.bond_c2 = np.delete(self.bond_c2, delete)

        print('Total new bonds: %d' % len(self.bond_c1))

    def make_initiators_real(self):

        # When c1 bonds to a c2 atom due to initiation, the c1 adjacent to c2 (not the c1 atom involved in bonding)
        # gains an initiator which is just represented as a hydrogen
        # If a radical c2 reacts, then there is no initiator that needs to be made real

        # number of atoms between c1 and c2 (based on order that they are written in the coordinate file)
        diff = np.array(self.xlink_residue.c2_atom_indices) - np.array(self.xlink_residue.c1_atom_indices)
        bond_c2 = np.array([x for x in self.bond_c2 if (x - 1) not in self.radicals])  # don't include radicals

        if len(bond_c2) > 0:
            tails = bond_c2 % len(self.xlink_residue.c2_atoms)  # define which tail each c2 atom is part of
            c1 = bond_c2 - diff[tails]  # indices of c1 atoms which gain a hydrogen atom

            # Get (serial) indices of hydrogen atoms to make real
            for i in c1:
                ndx = np.where(self.vsites[:, 1] == i)[0]  # index 1 contains serial carbon index used to construct vsite
                self.initiators.append(int(self.vsites[ndx, 0]))  # index 0 contains serial index of dummy H atom

            # add bonds -- this can be made more general
            for i, x in enumerate(c1):
                self.all_bonds.append([x, self.initiators[i]])
                self.xlink_residue_atoms.type[self.initiators[i] - 1] = 'hc'  # change from dummy hydrogen to real hydrogen
                self.xlink_residue_atoms.mass[self.initiators[i] - 1] = 1.008  # give mass to the former dummy atom
                self.xlink_residue_atoms.type[x - 1] = 'c3'  # make carbon sp3 hybridized

    def identify_new_radicals(self):

        # first remove reacted radicals
        self.radicals = [x for x in self.radicals if (x + 1) not in self.bond_c2]

        # c2 next to c1 that just bonded to a different c2 is now a radical
        if len(self.bond_c1) > 0:
            diff = np.array(self.xlink_residue.c1_atom_indices) - np.array(self.xlink_residue.c2_atom_indices)
            tails = self.bond_c1 % len(self.xlink_residue.c1_atoms)  # define which tail each c2 atom is part of
            for i in self.bond_c1 - diff[tails]:  # indices of c2 atoms which become radicals
                self.radicals.append(i - 1)
                self.xlink_residue_atoms.type[i - 1] = 'c2'  # change atom type from ce to c2

    def remove_improper_dihedrals(self):

        # keep improper dihedrals involving oxygen
        oxygen = self.xlink_residue.get_carbonyl_oxygens()

        rm_imp = []
        # all improper dihedrals that c1 is involved with do not involve the carbonyl oxygen
        for i in self.bond_c1:
            imps = np.where(self.impropers == i)[0]
            for j in imps:
                rm_imp.append(j)

        for i in self.bond_c2:
            imps = np.where(self.impropers == i)[0]
            for j in imps:
                test = [x for x in (self.impropers[j, :] % self.xlink_residue.natoms) if x in oxygen]
                if not test:  # if list is empty, then there is no carbonyl carbon, so we should remove the dihedral
                    rm_imp.append(j)

        rm_imp = np.unique(rm_imp)  # should be unnecessary

        self.impropers = np.delete(self.impropers, rm_imp, axis=0)

    def identify_terminated(self, rad_term_frac=0.5):

        if self.iteration > 1:  # don't do this on the first iteration

            # terminate enough radicals to keep the number of radicals stable
            nterm = int(rad_term_frac * len(self.radicals))
            chosen_terminated_radicals = np.random.choice(len(self.radicals), size=nterm, replace=False)
            self.terminated_radicals = np.array(self.radicals)[chosen_terminated_radicals]
            self.radicals = np.delete(self.radicals, chosen_terminated_radicals).tolist()

            # bond dummy hydrogen to terminated radicals
            for i in self.terminated_radicals:
                ndx = np.where(self.vsites[:, 1] == (i + 1))[0]  # index 1 contains serial carbon index used to construct vsite
                H = int(self.vsites[ndx, 0])  # index 0 contains serial index of dummy H atom
                self.all_bonds.append([i + 1, H])
                self.xlink_residue_atoms.type[H - 1] = 'hc'  # change from dummy hydrogen to real hydrogen
                self.xlink_residue_atoms.mass[H - 1] = 1.008  # give mass to the former dummy atom
                self.xlink_residue_atoms.type[i] = 'c3'  # make carbon sp3 hybridized
                self.initiators.append(H)  # add to initiators so that the virtual site gets deleted

            for i in self.terminated_radicals:
                self.terminate.append(i + 1)  # I think radicals are numbered from 0

        for i in self.bond_c1:
            self.terminate.append(i)
        for i in self.bond_c2:
            self.terminate.append(i)
            if i not in self.radicals:
                self.terminate.append(i + 1)  # terminate c1 adjacent to bonding c2 unless already terminated previously

    def remove_virtual_sites(self):

        # remove virtual sites for atoms that became real
        rm_vsite = []
        for i in range(self.vsites.shape[0]):
            if int(self.vsites[i, 0]) in self.initiators:
                rm_vsite.append(i)

        self.vsites = np.delete(self.vsites, rm_vsite, axis=0)  # delete vsite entries

    def bond(self):

        self.make_initiators_real()  # turn dummy atoms into real atoms
        self.identify_new_radicals()  # get indices of the latest radicals

        for i, x in enumerate(self.bond_c1):
            self.all_bonds.append([x, self.bond_c2[i]])
            # change atom types of newly bonded carbon atoms
            self.xlink_residue_atoms.type[x - 1] = 'c3'
            self.xlink_residue_atoms.type[self.bond_c2[i] - 1] = 'c3'

        self.nxlinks += len(self.bond_c1)  # update total number of cross-links in system

        self.identify_terminated(rad_term_frac=self.rad_frac_term)
        self.remove_virtual_sites()
        self.remove_improper_dihedrals()

        self.define_topology()  # make new angles, pairs and dihedrals lists

    def reload_coordinates(self, name):

        self.t = md.load(name)

    def update_lists(self):

        # reacted c1 and c2 atoms are no longer eligible for bonding
        # self.c1_atoms is numbered by index, while self.bond_c1 is number serially

        self.c1_atoms = [x for x in self.c1_atoms if (x + 1) not in self.terminate]
        self.c2_atoms = [x for x in self.c2_atoms if (x + 1) not in self.terminate]

        # re-initialize bonding carbon lists
        self.bond_c1 = []
        self.bond_c2 = []

        self.initiators = []

    def update_log(self):

        write_mode = 'a'
        if self.iteration == 1:
            write_mode = 'w'

        with open('xlink.log', '%s' % write_mode) as f:

            f.write('Iteration Number: %d \n' % self.iteration)
            f.write('Number of new cross-links: %d\n' % len(self.bond_c1))
            f.write('Total system cross-links: %d\n' % self.nxlinks)
            f.write('Cross-link density: %.2f %%\n' % (100. * (self.nxlinks / self.total_possible_xlinks)))
            f.write('+-----------------------------------+\n')
            f.write('| Serial indices of new cross-links |\n')
            f.write('+------+------+---------------------+\n')
            f.write('|{:^6}|{:^6}|\n'.format('c1', 'c2'))
            f.write('+------+------+\n')
            for i, x in enumerate(self.bond_c1):
                f.write('{:>5d} -- {:<5d}\n'.format(x, self.bond_c2[i]))
            f.write('Total radicals terminated: %.2d\n' % len(self.terminated_radicals))
            f.write('Total radicals left in system: %d\n' % len(self.radicals))
            f.write('Percent terminated: %.2f %%\n' % (100 * len(self.terminate) / (self.total_possible_terminated)))
            f.write('%s\n' % (80*'-'))

    def cleanup(self, out='xlinked.gro', dummy_atom_name='hc_d'):

        """
        Remove virtual sites in .gro and topology.
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


location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in

if __name__ == "__main__":

    args = initialize().parse_args()

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    start = time.time()
    initial_message = '# Cross-linking %s #' % args.initial
    print('#'*len(initial_message))
    print(initial_message)
    print('#'*len(initial_message))

    # if needed
    dt_short = 0.0005  # simulation time step (ps) (unstable if timestep is too large)
    length_short = .025  # (integer for now, so this is as low as it goes) length of short nvt simulation in picoseconds

    # get the system ready for cross-linking
    # Add dummy atoms create topology for entire assembly
    print('Initializing system and adding dummy atoms...', end='', flush=True)
    sys = System(args.initial, args.residue, args.dummy_residue, dummy_name=args.dummy_name,
                 radical_reaction_percentage=args.rad_percent, reaction_percentage=args.percent,
                 radical_termination_fraction=args.rad_frac_term)
    print('Done! %s created' % args.dummy_name)

    print('Generating input files: %s, %s, %s and %s ...' % (args.xlink_top_name, args.topname, args.mdp_em,
          args.mdp_nvt),  end='', flush=True)
    sys.write_assembly_topology()

    # generate input files that will be used throughtout cross-linking process
    full_top = SystemTopology(args.dummy_name, ff=args.forcefield, xlink=True, xlinked_top_name=args.xlink_top_name)  # topology with other residues and forcefield
    full_top.write_top(name=args.topname)
    mdpshort = SimulationMdp(args.dummy_name, T=args.temperature, em_steps=args.em_steps, time_step=args.dt,
        length=args.length, p_coupling='semiisotropic', barostat='berendsen', genvel='yes', restraints=False,
                             xlink=True, bcc=False)
    mdpshort.write_em_mdp(out='%s' % args.mdp_em.split('.')[0])
    mdpshort.write_nvt_mdp(out='%s' % args.mdp_nvt.split('.')[0])

    mdpshort = SimulationMdp(args.dummy_name, T=args.temperature, em_steps=args.em_steps, time_step=dt_short,
        length=length_short, p_coupling='semiisotropic', barostat='berendsen', genvel='yes', restraints=False,
                             xlink=True, bcc=False)
    mdpshort.write_nvt_mdp(out='nvt_short')
    print('Done!')

    # energy minimize starting configuration to get dummies in the right place
    print('Energy minimizing %s...' % args.dummy_name, end='', flush=True)
    sys.simulate(args.dummy_name, mdp=args.mdp_em, top=args.topname, out='em_%s' % args.dummy_name.split('.')[0],
                 parallel=args.parallelize, np=args.nproc)
    print('Done!')

    # Rest of iterations
    stagnated_iterations = 0
    while (len(sys.terminate) / sys.total_possible_terminated) < (args.density / 100.) and stagnated_iterations < args.stagnation:

        print('-'*80)
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

        # energy minimize then run short NVT simulation
        print('Energy minimizing new cross-links...', end='', flush=True)
        if sys.iteration == 1:
            sys.simulate('em_%s' % args.dummy_name.split('.')[0], mdp=args.mdp_em, top=args.topname, out='em',
                         parallel=args.parallelize, np=args.nproc)
        else:
            sys.simulate('nvt.gro', mdp=args.mdp_em, top=args.topname, out='em', parallel=args.parallelize,
                         np=args.nproc)
        print('Done!')

        print('Running %.1f ps NVT simulation...' % args.length, end='', flush=True)
        sys.simulate('em.gro', mdp=args.mdp_nvt, top=args.topname, out='nvt', parallel=args.parallelize, np=args.nproc)
        print('Done!')

        print('\nTotal new cross-links: %d' % len(sys.bond_c1))
        print('Total system cross-links: %d' % sys.nxlinks)
        print('Cross-link density: %.2f %%' % (100.*(sys.nxlinks / sys.total_possible_xlinks)))
        print('Total radicals terminated: %.2d' % len(sys.terminated_radicals))
        print('Total radicals left in system: %d' % len(sys.radicals))
        print('Percent terminated: %.2f' % (100 * len(sys.terminate) / sys.total_possible_terminated))

        if len(sys.bond_c1) > 0:
            stagnated_iterations = 0
        else:
            stagnated_iterations += 1

        # update log
        sys.update_log()
        sys.iteration += 1

    print('Removing remaining dummy atoms from .gro and .itp files...', end='', flush=True)
    sys.reload_coordinates('nvt.gro')
    sys.update_lists()
    sys.identify_terminated(rad_term_frac=1)  # terminate all radicals
    sys.remove_virtual_sites()
    sys.define_topology()
    sys.cleanup(out=args.output_gro)
    print('Done!')

    print('Energy minimizing %s...' % args.output_gro, end='', flush=True)
    sys.simulate(args.output_gro, mdp=args.mdp_em, top=args.topname, out=args.output_gro.split('.')[0],
                 parallel=args.parallelize, np=args.nproc)
    print('Done!')

    print('Finished cross-linking in %.2f seconds' % (time.time() - start))
