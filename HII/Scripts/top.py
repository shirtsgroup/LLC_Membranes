#! /usr/bin/env python

import Atom_props
import numpy as np


def get_indices(a, vsites):
    # find the indices of all fields that need to be modified
    atoms_index = 0  # find index where [ atoms ] section begins
    while a[atoms_index].count('[ atoms ]') == 0:
        atoms_index += 1

    bonds_index = 0  # find index where [ bonds ] section begins
    while a[bonds_index].count('[ bonds ]') == 0:
        bonds_index += 1

    pairs_index = 0  # find index where [ pairs ] section begins
    while a[pairs_index].count('[ pairs ]') == 0:
        pairs_index += 1
        if pairs_index == len(a):
            pairs_index = None
            break

    angles_index = 0  # find index where [ angles ] section begins
    while a[angles_index].count('[ angles ]') == 0:
        angles_index += 1
        if angles_index == len(a):
            angles_index = None
            break

    dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
    while a[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
        dihedrals_p_index += 1
        if dihedrals_p_index == len(a):
            dihedrals_p_index = None
            break

    dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
    while a[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
        dihedrals_imp_index += 1
        if dihedrals_imp_index == len(a):
            dihedrals_imp_index = None
            break

    if vsites:
        vsite_index = 0  # find index where [ dihedrals ] section begins (propers)
        while a[vsite_index].count('[ virtual_sites4 ]') == 0:
            vsite_index += 1
    else:
        vsite_index = None

        return atoms_index, bonds_index, pairs_index, angles_index, dihedrals_p_index, dihedrals_imp_index, vsite_index


class Top(object):

    def __init__(self, name, vsites=False):

        self.full = []
        with open(name, 'r') as f:
            for line in f:
                self.full.append(line)

        self.atoms_index, self.bonds_index, self.pairs_index, self.angles_index, self.dihedrals_p_index, \
        self.dihedrals_imp_index, self.vsite_index = get_indices(self.full, vsites)

        self.natoms = 0
        self.atoms = []
        self.atom_masses = []
        while self.full[self.atoms_index + self.natoms + 2] != '\n':
            self.atoms.append(str.strip(self.full[self.atoms_index + self.natoms + 2][5:17]))
            self.atom_masses.append(Atom_props.mass[str.strip(self.full[self.atoms_index + self.natoms + 2][5:17])])
            self.natoms += 1

        self.atom_masses = np.array(self.atom_masses)

        self.nbonds = 0
        while self.full[self.bonds_index + self.nbonds + 2] != '\n':
            self.nbonds += 1

    def find_indices(self, atoms):

        indices = []
        for i in range(self.atoms_index + 2, self.atoms_index + 2 + self.natoms):
            if str.strip(self.full[i][22:28]) in atoms:
                indices.append(i - self.atoms_index - 2)

        return indices

    def get_bonds_between(self, atoms1, atoms2):

        if type(atoms1) is not list: # allow single atoms to be passed without formatting as a list
            atoms1 = [atoms1]
        if type(atoms2) is not list:
            atoms2 = [atoms2]

        bonds = []

        for i in range(self.bonds_index + 2, self.bonds_index + 2 + self.nbonds):
            fields = [int(x) for x in self.full[i].split()]  # get all fields and convert into integers
            if fields[0] in atoms1:
                if fields[1] in atoms2:
                    bonds.append(fields[:2])
            if fields[0] in atoms2:
                if fields[1] in atoms1:
                    bonds.append(fields[:2])

        return bonds

    def get_bonds_to(self, atoms):
        """
        Get the indices of all the bonds to atoms
        :param atoms: Find all atoms bonded to atoms specified. (list of indices, or a single index)
        :return: dictionary where keys are atom indices from 'atoms' and values are bonds to that atom
        """

        if type(atoms) is not list:  # allow single atoms to be passed without formatting as a list
            atoms = [atoms]

        bonds = {}
        for i in atoms:
            bonds[i] = []  # initialize emtpy lists for each atom

        for i in range(self.bonds_index + 2, self.bonds_index + 2 + self.nbonds):
            fields = [int(x) for x in self.full[i].split()]  # get all fields and convert into integers
            if fields[0] in atoms:
                bonds[fields[0]].append(fields[1])
            if fields[1] in atoms:
                    bonds[fields[1]].append(fields[0])

        return bonds

    def get_type(self, atoms):
        """
        Get the atom types of 'atoms'
        :param atoms: indices of atoms whose atom types we want
        :return: atom types (list)
        """

        if type(atoms) is not list:
            atoms = [atoms]

        types = []
        for i in atoms:
            line = self.full[self.atoms_index + 1 + i]
            types.append(line.split()[1])

        return types

    def change_atom_type(self, atoms, new):
        """
        Change the atom types of the selected atoms to 'new'
        :param atoms: Indices of atoms whose types will be changed (list)
        :param new: Name of new atomtype (str)
        :return: modified topology
        """

        for i in atoms:
            ndx = i + self.atoms_index + 1  # plus 1 to convert between serial and index
            self.full[ndx] = self.full[ndx].replace(self.full[ndx][5:10], '{0: >5}'.format(new))

    def write_top(self, out):

        with open(out, 'w') as f:
            for line in self.full:
                f.write(line)

    def insert_atom(self, name, charge=0.0):

        ndx = self.atoms_index + self.natoms
        self.full.insert(ndx + 2, '{:5d}{:>5s}{:6d}{:>6s}{:>6s}{:>7d}{:13.6f}{:13.6f}\n'.format(ndx - self.atoms_index + 1, name,
        int(str.strip(self.full[ndx - 1][10:16])), str.strip(self.full[ndx - 1][16:22]), 'v',
        ndx - self.atoms_index + 1, charge, 0))

        self.natoms += 1

    def insert_vsite2(self, site, atom1, atom2, a):
        """
        Insert virtual site and/or virtual site section to exisiting topology
        :param site: index of virtual site atom
        :param atom1: 1st atom used to construct site
        :param atom2: 2nd atom used to construct site
        :param a: fraction of placement of vsite between atoms 1 and 2. Atom will be place a away from atom1 and 1-a
        from atom 2
        """

        if self.vsite_index == 0:
            self.full.append('[ virtual_sites2 ]\n')
            self.vsite_index = len(self.full)

        self.full.append('{:5d}{:6d}{:6d}{:6d}{:6.3f}\n'.format(site, atom1, atom2, 1, a))

    def insert_bond(self, atom1, atom2):

        ndx = self.bonds_index + self.nbonds + 2
        self.full.insert(ndx, '{:6d}{:7d}{:3d}\n'.format(atom1, atom2, 1))
        self.nbonds += 1