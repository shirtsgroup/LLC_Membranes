#!/usr/bin/env python

import os
import mdtraj as md
import Atom_props
import numpy as np
import sys

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

ions_mw = {}
ions_mw['NA'] = 22.99


class Residue(object):

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
            self.mw = ions_mw[name]
            self.mass = {}  # key = atom name, value = mass
            self.mass[name] = ions_mw[name]

        else:
            try:
                self.t = md.load('%s.pdb' % name,
                            standard_names=False)  # see if there is a solute configuration in this directory
            except OSError:
                try:
                    self.t = md.load('%s/../top/topologies/%s.pdb' % (script_location, name),
                                standard_names=False)  # check if the configuration is
                    # located with all of the other topologies
                except OSError:
                    raise OSError('No residue %s found. Perhaps you have not made a %s.pdb yet?' % (name, name))
                    #sys.exit('No residue %s found. Perhaps you have not made a %s.pdb yet?' % (name, name))

            try:
                f = open('%s.itp' % name, 'r')
            except FileNotFoundError:
                try:
                    f = open('%s/../top/topologies/%s.itp' % (script_location, name), 'r')
                except FileNotFoundError:
                    sys.exit('No topology %s.itp found' % name)

            self.res = [a.residue.name for a in self.t.topology.atoms]

            if len(set(self.res)) > 1:
                self.residues = []
                for i in set(self.res):
                    self.residues.append(Residue(i))

            else:
                self.resname = self.res[0]
                self.itp = []
                for line in f:
                    self.itp.append(line)

                f.close()

                atoms_index = 0
                while self.itp[atoms_index].count('[ atoms ]') == 0:
                    atoms_index += 1

                self.indices = {}  # key = atom name , value = index
                self.names = {}  # key = index, value = atom name
                self.mass = {}  # key = atom name, value = mass
                self.charges = {}

                atoms_index += 2
                self.natoms = 0
                while self.itp[self.natoms + atoms_index] != '\n':
                    data = self.itp[self.natoms + atoms_index].split()
                    self.indices[data[4]] = int(data[0])
                    self.names[int(data[0])] = data[4]
                    self.mass[data[4]] = float(data[7])
                    self.charges[data[4]] = float(data[6])
                    self.natoms += 1

                self.mw = sum(self.mass.values())

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

                self.bonds = {}
                for i in range(self.natoms):
                    self.bonds[i] = []
                    involvement = [x for x in bonds if i + 1 in x]
                    for pair in involvement:
                        atom = [x - 1 for x in pair if x != (i + 1)][0]
                        self.bonds[i].append(atom)


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
            self.mw = Atom_props.mass[name]
            self.charge = Atom_props.charge[name]
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
                self.mw += Atom_props.mass[a.name]

            self.com = np.zeros([3])  # center of mass of solute
            for i in range(self.xyz.shape[1]):
                self.com += self.xyz[0, i, :] * Atom_props.mass[self.names[i]]
            self.com /= self.mw


class LC(object):
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

        self.name = os.path.splitext(name)[0]  # build monomer name

        a = []
        with open('%s/../top/topologies/%s' % (script_location, name)) as f:
            for line in f:
                a.append(line)

        t = md.load("%s/../top/topologies/%s" % (script_location, name))
        self.LC_positions = t.xyz[0, :, :]
        self.LC_names = [a.name for a in t.topology.atoms]
        self.LC_residues = [a.residue.name for a in t.topology.atoms]

        self.full = a
        P = []
        P_ndx = []
        L = np.zeros([2], dtype=int)
        C1 = []
        C2 = []
        C1_ndx = []
        C2_ndx = []
        I_ndx = []
        carb = []
        PDA = []
        self.tail_atoms = []
        self.no_ions = 0
        self.ions = []
        self.MW = 0
        self.dummies = []

        if name.endswith('.gro'):

            # Lists of Plane atoms, line atoms, and dummy atoms
            residues = [str.strip(a[2][5:10])]
            nres = []
            res_count = 0
            benzene_carbons = []
            for i in range(2, len(a) - 1):
                res = str.strip(a[i][5:10])
                self.MW += Atom_props.mass[(str.strip(a[i][10:15]))]
                if res in residues:
                    res_count += 1
                else:
                    residues.append(res)
                    nres.append(res_count)
                    res_count = 1
                if a[i].count(';') != 0:
                    fields = a[i].split(';')
                    annotations = fields[1].split()
                    if 'P' in annotations:
                        P.append(str.strip(a[i][10:15]))
                        P_ndx.append(i - 2)
                    if 'L1' in annotations:
                        L[1] = i - 2
                    if 'L2' in annotations:
                        L[0] = i - 2
                    if 'R' in annotations:
                        self.ref_atom_index = i - 2
                    if 'C1' in annotations:
                        C1.append(str.strip(a[i][10:15]))
                        C1_ndx.append(int(a[i][15:20]))
                    if 'C2' in annotations:
                        C2.append(str.strip(a[i][10:15]))
                        C2_ndx.append(int(a[i][15:20]))
                    if 'I' in annotations:
                        self.no_ions += 1
                        ion = str.strip(a[i][10:15])
                        self.valence = Atom_props.charge[ion]
                        if ion not in self.ions:
                            self.ions.append(ion)
                        I_ndx.append(i - 2)
                    if 'B' in annotations:
                        benzene_carbons.append(str.strip(a[i][10:15]))
                    if 'C' in annotations:
                        carb.append(i - 2)
                    if 'PDA' in annotations:
                        PDA.append(str.strip(a[i][10:15]))
                    if 'T' in annotations:
                        self.tail_atoms.append(str.strip(a[i][10:15]))
                    if 'D' in annotations:
                        self.dummies.append(str.strip(a[i][10:15]))

            nres.append(res_count)

        elif name.endswith('.pdb'):

            # Lists of Plane atoms, line atoms, and dummy atoms
            start = 0
            while a[start].count('ATOM') == 0:
                start += 1
            end = start
            while a[end].count('ATOM') != 0:
                end += 1

            residues = [str.strip(a[start][17:22])]
            nres = []
            res_count = 0
            benzene_carbons = []
            for i in range(start, end):
                res = str.strip(a[i][17:22])
                if res in residues:
                    res_count += 1
                else:
                    residues.append(res)
                    nres.append(res_count)
                    res_count = 1
                if a[i].count(';') != 0:
                    fields = a[i].split(';')
                    annotations = fields[1].split()
                    if 'P' in annotations:
                        P.append(str.strip(a[i][11:16]))
                    if 'L' in annotations:
                        L.append(i - start)  # adjust for top lines and index (count from 0 rather than 1)
                    if 'R' in annotations:
                        self.ref_atom_index = i - start
                    if 'C1' in annotations:
                        C1.append(str.strip(a[i][11:16]))
                        C1_ndx.append(int(str.strip(a[i][0:5])))
                    if 'C2' in annotations:
                        C2.append(str.strip(a[i][11:16]))
                        C2_ndx.append(int(str.strip(a[i][0:5])))
                    if 'I' in annotations:
                        self.valence = Atom_props.charge[str.strip(a[i][11:16])]
                    if 'B' in annotations:
                        benzene_carbons.append(str.strip(a[i][11:16]))

            nres.append(res_count)

        self.natoms = sum(nres) - self.no_ions
        self.planeatoms = P
        self.plane_indices = P_ndx
        self.lineatoms = L
        self.residues = residues
        self.nresidues = nres
        self.c1_atoms = C1
        self.c2_atoms = C2
        self.c1_index = C1_ndx
        self.c2_index = C2_ndx
        self.benzene_carbons = benzene_carbons
        self.ion_indices = I_ndx
        self.carboxylate_indices = carb
        self.pore_defining_atoms = PDA

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
