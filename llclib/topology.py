#!/usr/bin/env python

import os
import mdtraj as md
import Atom_props
import numpy as np

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

        else:
            try:
                t = md.load('%s.pdb' % name,
                            standard_names=False)  # see if there is a solute configuration in this directory
            except OSError:
                try:
                    t = md.load('%s/../top/topologies/%s.pdb' % (script_location, name),
                                standard_names=False)  # check if the configuration is
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
                    exit()

            res = set([a.residue.name for a in t.topology.atoms])

            if len(set(res)) > 1:
                self.residues = []
                for i in set(res):
                    self.residues.append(Residue(i))

            else:

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
