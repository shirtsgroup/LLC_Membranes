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
from scipy.sparse import lil_matrix

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in


def initialize():

    parser = argparse.ArgumentParser(description='Crosslink LLC structure')  # allow input from user

    parser.add_argument('-i', '--initial', default='wiggle_init.gro', help = 'Name of input file')
    parser.add_argument('-b', '--build_mon', default='NAcarb11Vd', type=str, help='Monomer the system is built with')
    parser.add_argument('-l', '--layers', default=20, help='Number of layers')
    parser.add_argument('-p', '--pores', default=4, help='Number of pores')
    parser.add_argument('-n', '--no_monomers', default=5, help='Number of monomers per layer')
    parser.add_argument('-c', '--cutoff', default=5, type=float, help='Cutoff percentage for cross-linking. Bottom x %'
                                                              ' of the distribution of distances will be cross-linked ')
    parser.add_argument('-e', '--term_prob', default=5, type=float, help='Termination probability (%)')
    parser.add_argument('-y', '--topology', default='crosslinked_new.itp', help='Topology file to be analyzed and modified')
    parser.add_argument('-r', '--iteration', default=0, type=int, help='Iteration number of crosslinking process')
    parser.add_argument('-d', '--cutoff_rad', default=10, help='Cutoff Distance for radical reaction')
    parser.add_argument('-x', '--xlinks', default=0, help='Total number of c1-c2 crosslinks')
    parser.add_argument('-S', '--stop', default='no', help='Is crosslinking reaction finished')
    parser.add_argument('-m', '--monomer', default='monomer2', help='Which monomer was the structure built with')
    parser.add_argument('--nogap', help='System is prepared with no gap between membrane layers', action="store_true")

    return parser


class Residue(object):

    def __init__(self, name):

        self.name = name
        self.is_ion = False
        self.check_ion()
        self.virtual_sites = None
        self.c1_atoms = None
        self.c2_atoms = None
        self.c1_atom_indices = None
        self.c2_atom_indices = None

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

        with open('%s/../top/topologies/ions.txt' % location) as f:
            ions = []
            for line in f:
                if line[0] != '#':
                    ions.append(str.strip(line))

        if self.name in ions:
            self.is_ion = True

    def get_atoms(self):

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


class Topology():

    def __init__(self):

        self.residues = [Residue(a) for a in list(set([a.residue.name for a in self.t.topology.atoms]))]

        self.res_ndx = [i.name for i in self.residues].index(self.xlink_residue_name)
        self.xlink_residue = self.residues[self.res_ndx]
        self.xlink_residue.get_reactive_carbons()

        self.nresidues = {}
        for r in self.residues:
            self.nresidues[r.name] = self.all_residues.count(r.name) // r.natoms

        self.all_bonds = []
        natoms = self.residues[self.res_ndx].natoms
        bonds = self.residues[self.res_ndx].bonds
        for i in range(self.nresidues[self.xlink_residue_name]):
            for j in range(len(bonds)):
                self.all_bonds.append([int(bonds[j][0]) + i * natoms, int(bonds[j][1]) + i * natoms])

        # get pairs, angles, dihedrals
        self.pairs = None
        self.angles = None
        self.dihedrals = None
        self.adjacency_matrix = None

        self.define_topology()

        self.initiators = []  # indices of new initiator atoms which need to switched from virtual sites
        self.radicals = []  # indices of carbon atoms that are now radicals

    def define_topology(self):

        start = time.time()
        print("Creating adjacency matrix")
        self.create_adjacency_matrix()
        print("Finding angles")
        self.find_angles()
        print("Finding dihedrals")
        self.find_dihedrals()
        print("Getting Pairs")
        self.find_pairs()
        print('did that in %s seconds' % (time.time() - start))

    def write_assembly_topology(self, out='assembly.itp'):

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
                    f.write('{:>5d}{:>5s}{:>6s}{:>6s}{:>6s}{:>7s}{:>13s}{:>13s}\n'.format(int(atom_info[j][0]) + i*natoms,
                            atom_info[j][1], atom_info[j][2], atom_info[j][3], atom_info[j][4], atom_info[j][5],
                            atom_info[j][6], atom_info[j][7]))

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
            impropers = self.residues[self.res_ndx].improper_dihedrals

            for i in range(nres):
                for j in range(len(impropers)):
                    f.write('{:<6d}{:<7d}{:<7d}{:<7d}{:<7d}\n'.format(int(impropers[j][0]) + i * natoms,
                            int(impropers[j][1]) + i * natoms, int(impropers[j][2]) + i * natoms, int(impropers[j][3])
                            + i * natoms, 1))

            f.write('\n[ virtual_sites4 ]\n')
            f.write(';Site  from                         funct   a          b          d\n')
            vsites = self.residues[self.res_ndx].virtual_sites
            for i in range(nres):
                for j in range(len(vsites)):
                    if (int(vsites[j][0]) + i * natoms) not in self.initiators:
                        f.write('{:<8d}{:<6d}{:<6d}{:<6d}{:<8d}{:<8s}{:<11s}{:<11s}{:<11s}\n'.format(int(vsites[j][0]) +
                                i * natoms, int(vsites[j][1]) + i * natoms, int(vsites[j][2]) + i * natoms,
                                int(vsites[j][3]) + i * natoms, int(vsites[j][4]) + i * natoms, vsites[j][5], vsites[j][6],
                                vsites[j][7], vsites[j][8]))

    def create_adjacency_matrix(self):

        natoms = self.nresidues[self.xlink_residue_name] * self.residues[self.res_ndx].natoms
        #self.adjacency_matrix = np.zeros([natoms, natoms], dtype=bool)  # 27 --> 3 GB when changed dtype to bool
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
                    self.angles.append((i, bond[0], bond[1]))

    def find_dihedrals(self):

        self.dihedrals = []
        for angle in self.angles:
            # nonzero = np.nonzero(self.adjacency_matrix[angle[0] - 1])[0] + 1
            #nonzero = self.adjacency_matrix[angle[0] - 1].nonzero()[1] + 1
            nonzero = np.array(self.adjacency_matrix[angle[0] - 1].rows[0]) + 1
            for i in nonzero:
                if i < angle[2] and i not in angle:
                    self.dihedrals.append((i, angle[0], angle[1], angle[2]))

    def find_pairs(self):

        self.pairs = []
        for i in self.dihedrals:
            self.pairs.append((i[0], i[3]))


class System(Topology):

    def __init__(self, initial_configuration, residue, dummy_residue, dummy_name='dummies.gro'):

        add_dummies(md.load(initial_configuration), residue, dummy_residue,
                    out=dummy_name)  # add dummy atoms to the initial configuration

        self.xlink_residue_name = dummy_residue
        self.t = md.load(dummy_name)
        self.all_residues = [a.residue.name for a in self.t.topology.atoms]

        super().__init__()  # inherit from Topology()

        self.c1_atoms = [a.index for a in self.t.topology.atoms if a.name in self.xlink_residue.c1_atoms]
        self.c2_atoms = [a.index for a in self.t.topology.atoms if a.name in self.xlink_residue.c2_atoms]

        self.bond_c1 = None
        self.bond_c2 = None

    def energy_minimize(self, configuration, mdp='em.mdp', top='topol.top', out='em'):
        """
        Energy minimize a configuration using existing .mdp files
        """

        p1 = subprocess.Popen(
            ["gmx", "grompp", "-p", "%s" % top, "-f", "%s" % mdp, "-o", "%s" %out, "-c", "%s" % configuration],
            stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # generate atomic level input file
        p1.wait()

        p2 = subprocess.Popen(["gmx", "mdrun", "-deffnm", "%s" % out], stdout=open(os.devnull, 'w'),
                              stderr=subprocess.STDOUT)  # run energy minimization
        p2.wait()

    def select_eligible_carbons(self, percent=1, cut=0.6):
        """
        - Calculate the distance between atoms of interest. Find out which pairs meet the distance cutoff criteria
        - Head-to-tail addition dominates due to sterics and radical stabilization
        :param percent: the bottom percentage of the distribution of distances (below cut) that will be chosen for cross-linking
        :param cut: any distance beyond 'cut' will not be accounted for when calculating percent
        """
        # calculate minimum image distance between c1 and c2 atoms
        dim = self.t.unitcell_lengths[0]

        # shift so box is cubic temporarily. The minimum image distances are just an approximation using this approach
        theta = np.pi / 3
        self.t.xyz[..., 1] /= np.sin(theta)
        self.t.xyz[..., 0] = self.t.xyz[..., 0] - self.t.xyz[..., 1]*np.cos(theta)

        c1_len = len(self.c1_atoms)
        c2_len = len(self.c2_atoms)

        diff = np.zeros([c1_len*c2_len, 3])
        for i in range(c1_len):
            diff[i*c2_len: (i + 1)*c2_len] = self.t.xyz[0, self.c1_atoms[i], :] - self.t.xyz[0, self.c2_atoms, :]

        # make sure we are looking at minimum image distance
        diff = np.where(diff > 0.5 * dim, diff - dim, diff)
        diff = np.where(diff < -0.5 * dim, diff + dim, diff)

        # calculate distance between all relevant carbon atoms
        distance = np.linalg.norm(diff, axis=1)
        diff_ndx = np.argsort(distance)  # index of distances sorted from low to high
        eligible_c1 = (diff_ndx // c1_len)[c1_len:]  # exclude c1 and c2 atoms that are bonded to each other
        eligible_c2 = (diff_ndx % c1_len)[c2_len:]

        nc1 = len(self.xlink_residue.c1_atoms)
        nc2 = len(self.xlink_residue.c2_atoms)

        eligible_c1 = (eligible_c1 // nc1) * self.xlink_residue.natoms + \
                            np.array(self.xlink_residue.c1_atom_indices)[(eligible_c1 % nc1)]
        eligible_c2 = (eligible_c2 // nc2) * self.xlink_residue.natoms + \
                            np.array(self.xlink_residue.c2_atom_indices)[(eligible_c2 % nc2)]

        d = np.linalg.norm(self.t.xyz[-1, eligible_c1 - 1, :] - self.t.xyz[-1, eligible_c2 - 1, :], axis=1)
        cut_ndx = np.argmin(np.abs(d - cut))  # index where cut is exceeded

        # The final list of carbon atoms to be bonded
        self.bond_c1 = eligible_c1[:(int(percent * cut_ndx / 100))]
        self.bond_c2 = eligible_c2[:(int(percent * cut_ndx / 100))]

    def make_initiators_real(self):

        # When c1 bonds to a c2 atom due to initiation, the c1 adjacent to c2 (not the c1 atom involved in bonding)
        # gains an initiator which is just represented as a hydrogen

        # number of atoms between c1 and c2 (based on order that they are written in the coordinate file)
        diff = np.array(self.xlink_residue.c2_atom_indices) - np.array(self.xlink_residue.c1_atom_indices)
        tails = self.bond_c2 % len(self.xlink_residue.c2_atoms)  # define which tail each c2 atom is part of
        c1 = self.bond_c2 - diff[tails]  # indices of c1 atoms which gain a hydrogen atom

        # This loop can probably be made better if I set up the virtual site list differently
        # Get indices of hydrogen atoms to make real
        for i, ndx in enumerate(c1 % self.xlink_residue.natoms):
            j = 0
            # the following assumes that index 1 contains the index of the atom to which hydrogen will be bonded
            # It is the index of the atom bonded to the virtual site hydrogen
            while int(self.xlink_residue.virtual_sites[j][1]) != ndx:
                j += 1

            H_serial = int(self.xlink_residue.virtual_sites[j][0]) + self.xlink_residue.natoms * \
                       (c1[i] // self.xlink_residue.natoms)
            self.initiators.append(H_serial)

            # TODO: This doesn't work since topology only describes a single monomer. So this will need to be written
            # into the topology writer based on all of the new hydrogens
            # self.xlink_residue.atom_info[H_serial - 1][1] = 'hc'  # change atom type of hydrogen
            # self.xlink_residue.atom_info[H_serial - 1][-1] = 1.008  # given atom mass

        # add bonds
        for i, x in enumerate(c1):
            self.all_bonds.append([x, self.initiators[i]])

    def identify_new_radicals(self):

        # c2 next to c1 that just bonded to a different c2 is now a radical
        diff = np.array(self.xlink_residue.c1_atom_indices) - np.array(self.xlink_residue.c2_atom_indices)
        tails = self.bond_c1 % len(self.xlink_residue.c1_atoms)  # define which tail each c2 atom is part of
        for i in self.bond_c1 - diff[tails]:  # indices of c1 atoms which gain a hydrogen atom
            self.radicals.append(i)

    def bond(self):

        self.make_initiators_real()  # turn dummy atoms into real atoms
        self.identify_new_radicals()  # get indices of the latest radicals

        for i, x in enumerate(self.bond_c1):
            self.all_bonds.append([x, self.bond_c2[i]])

        self.define_topology()  # make new angles, pairs and dihedrals lists


if __name__ == "__main__":

    args = initialize().parse_args()

    # TODO: put the following into initialize() so they can be parsed by argparse
    dummy_name = 'dummies.gro'  # name of .gro file after dummy atoms are added
    xlink_top_name = 'assembly.itp'  # name of .itp file used to describe xlinked assembly
    topname = 'topol.top'  # name of full topology
    T = 300  # temperature to run system at
    em_steps = 5000  # number of steps for energy minimization
    dt = 0.001  # simulation time step
    length = 10  # length of simulation in picoseconds
    forcefield = 'gaff'
    mdp_em = 'em.mdp'
    mdp_nvt = 'nvt.mdp'

    # get the system ready for cross-linking
    # TODO: reorder this so that dummy configuration is energy minimized before any cross-linking
    sys = System(args.initial, 'HII', 'HIId', dummy_name=dummy_name)
    sys.select_eligible_carbons(percent=args.cutoff)
    sys.bond()
    sys.write_assembly_topology()
    # exit()
    # sys.write_assembly_topology()
    full_top = SystemTopology(dummy_name, ff='gaff', xlink=True)  # topology with other residues and forcefield
    full_top.write_top(name=topname, crosslinked_top_name=xlink_top_name)
    # generate input files that will be used throughtout process
    mdp = SimulationMdp(dummy_name, T=T, em_steps=em_steps, time_step=dt, length=length, p_coupling='semiisotropic',
                        barostat='berendsen', genvel='yes', restraints=False, xlink=True, bcc=False)
    mdp.write_em_mdp(out='%s' % mdp_em.split('.')[0])
    mdp.write_nvt_mdp(out='%s' % mdp_nvt.split('.')[0])

    # energy minimize starting configuration
    sys.energy_minimize(dummy_name, mdp=mdp_em, top=topname, out='em')

    """
    TODO: So far, this script will read in an initial configuration, add dummy atoms, create a topology by generating
    dihedrals and then energy minimize it. I'm in the middle of working on calculating the distance between reactive
    carbons. 
    
    After that, I'll need to create create new bonds, change atom-types where needed -- could probably build atom-types
    into topology rather than replacing strings.
    
    I'll also implement radicals.
    """
