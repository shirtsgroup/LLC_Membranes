#!/usr/bin/env python

import os
import mdtraj as md

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


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

            itp = []
            for line in f:
                itp.append(line)

            f.close()

            self.natoms = t.n_atoms

            atoms_index = 0
            while itp[atoms_index].count('[ atoms ]') == 0:
                atoms_index += 1

            self.indices = {}  # key = atom name , value = index
            self.names = {}  # key = index, value = atom name
            self.mass = {}  # key = atom name, value = mass

            atoms_index += 2
            for i in range(self.natoms):
                data = itp[i + atoms_index].split()
                self.indices[data[4]] = int(data[0])
                self.names[int(data[0])] = data[4]
                self.mass[data[4]] = float(data[7])

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

            self.bonds = {}
            for i in range(self.natoms):
                self.bonds[i] = []
                involvement = [x for x in bonds if i + 1 in x]
                for pair in involvement:
                    atom = [x - 1 for x in pair if x != (i + 1)][0]
                    self.bonds[i].append(atom)