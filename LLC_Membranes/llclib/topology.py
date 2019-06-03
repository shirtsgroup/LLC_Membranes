#!/usr/bin/env python

import os
import mdtraj as md
from LLC_Membranes.llclib import atom_props
import numpy as np
import sys
import subprocess

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

ions_mw = {}
ions_mw['NA'] = 22.99


class ReadItp(object):
    """ Read a GROMACS topology file and extract desired information

    :param name: name of GROMACS topology file

    :return: charge,
    """

    def __init__(self, name):

        try:
            f = open('%s.itp' % name, 'r')
        except FileNotFoundError:
            try:
                f = open('%s/../top/topologies/%s.itp' % (script_location, name), 'r')
            except FileNotFoundError:
                raise FileNotFoundError('No topology %s.itp found' % name)

        self.itp = []
        for line in f:
            self.itp.append(line)

        f.close()

        # Intialize atom properties
        self.natoms = 0
        self.indices = {}  # key = atom name , value = index
        self.names = {}  # key = index, value = atom name
        self.mass = {}  # key = atom name, value = mass
        self.charges = {}

        # Initialize annotation lists
        self.hbond_H = []  # hydrogen atoms capable of hbonding
        self.hbond_D = []  # hydrogen bond donor atoms
        self.hbond_A = []  # hydrogen bond acceptors
        self.residues = []
        self.planeatoms = []
        self.plane_indices = []
        self.benzene_carbons = []
        #self.lineatoms = np.zeros([2], dtype=int)
        self.lineatoms = [[], []]
        self.ref_atom_index = []
        self.c1_atoms = []
        self.c2_atoms = []
        self.c1_index = []
        self.c2_index = []
        self.ion_indices = []
        self.tail_atoms = []
        self.no_ions = 0
        self.ions = []
        self.MW = 0
        self.dummies = []
        self.valence = 0
        self.carboxylate_indices = []
        self.pore_defining_atoms = []

        # bonds
        self.bonds = {}

    def atoms(self, annotations=False):

        atoms_index = 0
        while self.itp[atoms_index].count('[ atoms ]') == 0:
            atoms_index += 1

        atoms_index += 2

        while self.itp[self.natoms + atoms_index] != '\n':

            data = self.itp[self.natoms + atoms_index].split()

            _, type, _, resname, atom_name, _, _, _ = data[:8]

            ndx = int(data[0]) - 1
            charge = float(data[6])

            try:
                mass = float(data[7])
            except ValueError:  # Throws error if annotation is present with semi colon directly next to mass
                mass = float(data[7].split(';')[0])

            self.indices[atom_name] = ndx
            self.names[ndx] = atom_name
            self.mass[atom_name] = float(mass)
            self.charges[atom_name] = float(charge)
            self.MW += mass

            if resname not in self.residues:
                self.residues.append(resname)

            if annotations:

                try:  # check for annotations on atom

                    annotations = self.itp[self.natoms + atoms_index].split(';')[1].split()

                    if 'H' in annotations:
                        self.hbond_H.append(atom_name)
                    if 'D' in annotations:
                        self.hbond_D.append(atom_name)
                    if 'A' in annotations:
                        self.hbond_A.append(atom_name)
                    if 'P' in annotations:
                        self.planeatoms.append(atom_name)
                        self.plane_indices.append(ndx)
                    if 'L1' in annotations:
                        self.lineatoms[1].append(ndx)
                    if 'L2' in annotations:
                        self.lineatoms[0].append(ndx)
                    if 'R' in annotations:
                        self.ref_atom_index.append(ndx)
                    if 'C1' in annotations:
                        self.c1_atoms.append(atom_name)
                        self.c1_index.append(ndx)
                    if 'C2' in annotations:
                        self.c2_atoms.append(atom_name)
                        self.c2_index.append(ndx)
                    if 'I' in annotations:
                        self.no_ions += 1
                        self.valence = charge
                        if atom_name not in self.ions:
                            self.ions.append(atom_name)
                        self.ion_indices.append(ndx)
                    if 'B' in annotations:
                        self.benzene_carbons.append(atom_name)
                    if 'C' in annotations:
                        self.carboxylate_indices.append(ndx)
                    if 'PDA' in annotations:
                        self.pore_defining_atoms.append(atom_name)
                    if 'T' in annotations:
                        self.tail_atoms.append(atom_name)
                    if 'D' in annotations:
                        self.dummies.append(atom_name)
                except IndexError:
                    pass

            self.natoms += 1

        #self.lineatoms = [np.array(x) for x in self.lineatoms]

    def organize_bonds(self):

        # find the bonds section
        bonds_index = 0
        while self.itp[bonds_index].count('[ bonds ]') == 0:
            bonds_index += 1
        bonds_index += 2

        bonds = []
        while self.itp[bonds_index] != '\n':
            bond_data = str.split(self.itp[bonds_index])[:2]
            bonds.append([int(bond_data[0]), int(bond_data[1])])
            bonds_index += 1

        for i in range(self.natoms):
            self.bonds[i] = []
            involvement = [x for x in bonds if i + 1 in x]
            for pair in involvement:
                atom = [x - 1 for x in pair if x != (i + 1)][0]
                self.bonds[i].append(atom)


class Residue(ReadItp):

    def __init__(self, name):

        self.name = name

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
            self.MW = ions_mw[name]
            self.mass = {}  # key = atom name, value = mass
            self.mass[name] = ions_mw[name]

        else:

            super().__init__(name)
            self.atoms(annotations=True)
            self.organize_bonds()
            # try:
            #     self.t = md.load('%s.pdb' % name,
            #                      standard_names=False)  # see if there is a solute configuration in this directory
            # except OSError:
            #     try:
            #         self.t = md.load('%s/../top/topologies/%s.pdb' % (script_location, name),
            #                          standard_names=False)  # check if the config is with all of the other topologies
            #     except OSError:
            #         raise OSError('No residue %s found. Perhaps you have not made a %s.pdb yet?' % (name, name))
            #
            # try:
            #     f = open('%s.itp' % name, 'r')
            # except FileNotFoundError:
            #     try:
            #         f = open('%s/../top/topologies/%s.itp' % (script_location, name), 'r')
            #     except FileNotFoundError:
            #         raise FileNotFoundError('No topology %s.itp found' % name)
            #
            # self.res = [a.residue.name for a in self.t.topology.atoms]
            #
            # if len(set(self.res)) > 1:
            #     self.residues = []
            #     for i in set(self.res):
            #         self.residues.append(Residue(i))
            #
            # else:
            #     self.resname = self.res[0]
            #     self.itp = []
            #     for line in f:
            #         self.itp.append(line)
            #
            #     f.close()
            #
            #     atoms_index = 0
            #     while self.itp[atoms_index].count('[ atoms ]') == 0:
            #         atoms_index += 1
            #
            #     self.indices = {}  # key = atom name , value = index
            #     self.names = {}  # key = index, value = atom name
            #     self.mass = {}  # key = atom name, value = mass
            #     self.charges = {}
            #
            #     self.hbond_H = []  # hydrogen atoms capable of hbonding
            #     self.hbond_D = []  # hydrogen bond donor atoms
            #     self.hbond_A = []  # hydrogen bond acceptors
            #
            #     atoms_index += 2
            #     self.natoms = 0
            #     while self.itp[self.natoms + atoms_index] != '\n':
            #
            #         data = self.itp[self.natoms + atoms_index].split()
            #         self.indices[data[4]] = int(data[0])
            #         self.names[int(data[0])] = data[4]
            #         self.mass[data[4]] = float(data[7])
            #         self.charges[data[4]] = float(data[6])
            #
            #         try:  # check for annotations on atom
            #             annotations = self.itp[self.natoms + atoms_index].split(';')[1].split()
            #             if 'H' in annotations:
            #                 self.hbond_H.append(data[4])
            #             if 'D' in annotations:
            #                 self.hbond_D.append(data[4])
            #             if 'A' in annotations:
            #                 self.hbond_A.append(data[4])
            #         except IndexError:
            #             pass
            #
            #         self.natoms += 1
            #
            #     self.mw = sum(self.mass.values())

                # # find the bonds section
                # bonds_index = 0
                # while self.itp[bonds_index].count('[ bonds ]') == 0:
                #     bonds_index += 1
                # bonds_index += 2
                #
                # bonds = []
                # while self.itp[bonds_index] != '\n':
                #     bond_data = str.split(self.itp[bonds_index])[:2]
                #     bonds.append([int(bond_data[0]), int(bond_data[1])])
                #     bonds_index += 1
                #
                # self.bonds = {}
                # for i in range(self.natoms):
                #     self.bonds[i] = []
                #     involvement = [x for x in bonds if i + 1 in x]
                #     for pair in involvement:
                #         atom = [x - 1 for x in pair if x != (i + 1)][0]
                #         self.bonds[i].append(atom)


class Molecule(object):

    def __init__(self, name):

        self.is_ion = False
        # check if residue is an ion
        with open('%s/../top/topologies/ions.txt' % script_location) as f:
            ions = []
            for line in f:
                if line[0] != '#':
                    ions.append(str.strip(line))

        if name in ions:
            self.is_ion = True
            self.residues = [name]
            self.names = [name]
            self.xyz = np.zeros([1, 1, 3])
            self.xyz[0, 0, :] = [0, 0, 0]
            self.natoms = 1
            self.mw = atom_props.mass[name]
            self.charge = atom_props.charge[name]
            self.resname = name
        else:
            try:
                t = md.load('%s.pdb' % name, standard_names=False)  # see if there is a solute configuration in this directory
            except OSError:
                try:
                    t = md.load('%s/../top/topologies/%s.pdb' % (script_location, name), standard_names=False)  # check if the configuration is
                    # located with all of the other topologies
                except OSError:
                    print('No residue %s found' % name)
                    exit()

            try:
                f = open('%s.itp' % name, 'r')
            except FileNotFoundError:
                try:
                    f = open('%s/../top/topologies/%s.itp' % (script_location, name), 'r')
                except FileNotFoundError:
                    print('No topology %s.itp found' % name)

            itp = []
            for line in f:
                itp.append(line)

            f.close()

            self.natoms = t.n_atoms

            atoms_index = 0
            while itp[atoms_index].count('[ atoms ]') == 0:
                atoms_index += 1

            atoms_index += 2
            self.charge = 0
            for i in range(self.natoms):
                self.charge += float(itp[atoms_index + i].split()[6])

            self.residues = [a.residue.name for a in t.topology.atoms]
            self.resname = self.residues[0]
            self.names = [a.name for a in t.topology.atoms]
            self.xyz = t.xyz

            self.mw = 0  # molecular weight (grams)
            for a in t.topology.atoms:
                self.mw += atom_props.mass[a.name]

            self.com = np.zeros([3])  # center of mass of solute
            for i in range(self.xyz.shape[1]):
                self.com += self.xyz[0, i, :] * atom_props.mass[self.names[i]]
            self.com /= self.mw


class LC(ReadItp):
    """A Liquid Crystal monomer has the following attributes which are relevant to building and crosslinking:

    Attributes:

        Description of annotations:
        "R" : reference atom: This atom defines the pore radius, r. It will be placed r nm from pore center
        "P" : plane atoms: 3 atoms defining a plane within the monomer which you want to be parallel to the xy plane
        "L" : line atoms: 2 atoms used to rotate monomers on xy plane so that the line created by line atoms goes
        through the pore center.
        "C1" : terminal vinyl carbon on tails. (for cross-linking)
        "C2" : second to last vinyl carbon on tails (for cross-linking)
        "B" : carbon atoms making up benzene ring

        name: A string representing the monomer's name.
        natoms: An integer accounting for the number of atoms in a single monomer.
        build_mon: Monomer used to build the unit cell
        images: Number of periodic images to be used in calculations
        c1_atoms: A list of atoms which will be involved in crosslinking as 'c1' -- See xlink.py
        c2_atoms: A list of atoms which will be involved in crosslinking as 'c2' -- See xlink.py
        tails: Number of tails on each monomer
        residues: A list of the minimum residue names present in a typical structure
        no_vsites: A string indicating whether there are dummy atoms associated with this monomer.

    Notes:
        Name of .gro and .itp are assumed to be the same unless otherwise specified. Whatever you pass to this class
        should be the name of the .gro/.itp file and it will read the annotations and directives
    """

    def __init__(self, name):

        super().__init__(name)

        self.atoms(annotations=True)

        self.name = name

        a = []
        with open('%s/../top/topologies/%s.gro' % (script_location, name)) as f:
            for line in f:
                a.append(line)

        t = md.load("%s/../top/topologies/%s.gro" % (script_location, name))
        self.LC_positions = t.xyz[0, :, :]
        self.LC_names = [a.name for a in t.topology.atoms]
        self.LC_residues = [a.residue.name for a in t.topology.atoms]

        # Things ReadItp gets wrong because it doesn't include the ion .itps
        self.natoms = len(self.LC_names)
        self.residues = np.unique(self.LC_residues)

        self.full = a

    def get_index(self, name):
        """
        Name of atoms whose index you want
        :param name: name listed in .gro file in 3rd column
        :return: index (serial) of the atom you want
        """
        ndx = -2
        for i in self.full:
            ndx += 1
            if str.strip(i[10:15]) == name:
                break

        return ndx


class Solute(Residue):

    def __init__(self, name):

        super().__init__(name)

        self.gro = []
        with open('%s/../top/topologies/%s.gro' % (script_location, name), 'r') as f:
            for line in f:
                self.gro.append(line)

        # vector defining the direction of solute. First entry is back of vector, second is front of vector
        self.direction_vector = [[], []]

        for i in range(2, len(self.gro) - 1):
            if self.gro[i].count(';') > 0:
                annotations = self.gro[i].split(';')[1].split()
                if annotations.count('Vback') > 0:
                    self.direction_vector[0].append(str.strip(self.gro[i][10:15]))
                if annotations.count('Vfront') > 0:
                    self.direction_vector[1].append(str.strip(self.gro[i][10:15]))


def map_atoms(indices, nres_atoms=1):
    """ Map the indices of a sub-system to indices of the full system

    :param indices: indices of atoms to map with respect to full system
    :param nres_atoms: number of atoms per residue

    :type indices: list
    :type nres_atoms: int

    :return: dictionary of mapped indices
    """

    index_map = {}
    nres = len(indices) // nres_atoms

    for i in range(nres):
        index_map[i] = indices[i*nres_atoms:(i + 1)*nres_atoms]

    return index_map


def fix_names(gro):
    """ Workaround for mdtraj. Fix atom names so they are the same as those shown in the gro file.

    :param gro: name of .gro file with all atom names in it

    :type gro: str
    """

    if gro.endswith('.gro'):
        gro = gro.split('.')[0]

    # if .pdb already exists, don't both remaking it -- this could be dangerous
    if not os.path.isfile('%s.pdb' % gro):

        convert_to_pdb = "gmx editconf -f %s.gro -o %s.pdb" % (gro, gro)
        p = subprocess.Popen(convert_to_pdb.split(), stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        p.wait()

    t = md.load('%s.pdb' % gro, standard_names=False)  # load doesn't have standard_names functionality for trr or xtc

    return [a.name for a in t.topology.atoms]
