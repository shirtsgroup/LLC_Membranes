#!/usr/bin/env python

import mdtraj as md
import numpy as np
import os
from LLC_Membranes.analysis import Atom_props

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


class Residue(object):

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
                # see if there is a solute configuration in this directory
                t = md.load('%s.pdb' % name, standard_names=False)
            except OSError:
                try:
                    # see if there is a solute configuration in this directory located with all of the other topologies
                    t = md.load('%s/../top/topologies/%s.pdb' % (script_location, name), standard_names=False)
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
            self.masses = []
            for i in range(self.xyz.shape[1]):
                self.masses.append(Atom_props.mass[self.names[i]])
                self.com += self.xyz[0, i, :] * self.masses[i]
            self.com /= self.mw
