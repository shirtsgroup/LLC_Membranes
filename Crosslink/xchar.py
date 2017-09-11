#!/usr/bin/env python

# Characterize Crosslinked systems

from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import argparse
import math
import os
import mdtraj as md
import lc_class
import numpy as np

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def initialize():
    parser = argparse.ArgumentParser(description = 'Determine where the crosslinks are')  # allow input from user

    # Flags
    parser.add_argument('-g', '--gro', default='wiggle_no_dummies.gro', help='Name of the final, cleaned .gro file after \
                                                                                crosslinking finishes')
    parser.add_argument('-t', '--top', default='crosslinked.itp', help='.itp file associated with input .gro file')
    parser.add_argument('-m', '--monomer', default='NAcarb11Vd.gro', help='Type of monomer used in crosslinked system')
    parser.add_argument('-p', '--no_pores', default=4, type=int, help='Number of Pores')
    parser.add_argument('-ntails', default=3, type=int, help='Number of tails on each monomer')
    parser.add_argument('-mpl', default=6, type=int, help='Monomers per layer')
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of layers')

    args = parser.parse_args()

    return args


class Xlinks:

    def __init__(self, monomer, structure):
        """
        Learn system details
        :param monomer: Name of monomer used to create full system
        :param structure: Name of crosslinked structure we are analysing
        """
        self.name = structure  # name of full system

        lc = lc_class.LC(monomer)
        self.c1_atoms = lc.c1_atoms  # name of c1 atoms
        self.c2_atoms = lc.c2_atoms  # name of c2 atoms
        self.c1_ndx = lc.c1_index  # indices of c1 atoms
        self.c2_ndx = lc.c2_index  # indices of c2 atoms
        self.mon_atoms = lc.natoms  # number of atoms in monomer
        self.c_spacing = np.array(self.c1_ndx) - np.array(self.c2_ndx)  # spacing of c atom indices

        t = md.load(structure)

        self.natoms = t.n_atoms  # number of atoms in full assembly
        self.c1 = [a.index + 1 for a in t.topology.atoms if a.name in self.c1_atoms]  # serial number of all c1 atoms
        self.c2 = [a.index + 1 for a in t.topology.atoms if a.name in self.c2_atoms]  # serial number of all c2 atoms
        self.pos = t.xyz  # positions of all atoms

        self.xlinks = []  # a list to hold all of the crosslinked pairs

        # since the number of atoms in the system has changed, each monomer may have slightly different number of atoms
        # To see which crosslinks are between atoms in the same pore, same monomer, same layer etc. we reindex them
        # so that only carbons are numbered from 0 to ntails*2*nmonomers
        self.reindexed = {}

        count = 0
        for i in range(self.natoms):
            if i in self.c1 or i in self.c2:
                self.reindexed[i] = count
                count += 1

    def is_xlink(self, atom, bonds):

        for i in bonds:  # check every atom that 'atom' is bonded to
            if atom in self.c1:  # if 'atom' is a c1 then check if it is bonded to any c2 atom from a different tail
                if i in self.c2 and (atom - i) not in self.c_spacing:  # should be more specific here. Need to figure out what kind of c2 it is (i.e. C34 or C17 etc)
                    if [i, atom] not in self.xlinks:
                        self.xlinks.append([atom, i])
            if atom in self.c2:
                if i in self.c1 and (i - atom) not in self.c_spacing:
                    if [i, atom] not in self.xlinks:
                        self.xlinks.append([atom, i])


class Topology:

    def __init__(self, top, system):

        self.name = top

        with open(top, 'r') as itp:

            a = []
            for l in itp:
                a.append(l)

        bonds_index = 0
        while a[bonds_index].count('[ bonds ]') == 0:
            bonds_index += 1

        bonds = [[] for _ in range(system.natoms)]

        bonds_index += 2
        while a[bonds_index] != '\n':
            b1 = int(a[bonds_index][0:6])
            b2 = int(a[bonds_index][6:13])
            bonds[b1 - 1].append(b2)
            bonds[b2 - 1].append(b1)
            bonds_index += 1

        self.bonds = bonds

if __name__ == "__main__":

    args = initialize()

    system = Xlinks(args.monomer, args.gro)  # get system information
    top = Topology(args.top, system)  # get topology information - bonds

    for i in range(len(top.bonds)):
        system.is_xlink(i + 1, top.bonds[i])  # identify all the cross links

    nxlinks = len(system.xlinks)

    # find out characteristics about the crosslinks - i.e. what relationship the xlinked carbons have with each other
    intra_mon = 0  # same monomer
    intra_layer = 0  # same layer
    intra_pore = 0  # same pore
    inter_pore = 0  # other pores
    for i in range(nxlinks):
        a = system.reindexed[system.xlinks[i][0]]
        b = system.reindexed[system.xlinks[i][1]]
        # [monomer no., layer no., pore no.]
        atom1 = [int(a / (2*args.ntails)), int(a / (2*args.ntails*args.mpl)), int(a/ (2*args.ntails*args.mpl*args.layers))]
        atom2 = [int(b / (2*args.ntails)), int(b / (2*args.ntails*args.mpl)), int(b/ (2*args.ntails*args.mpl*args.layers))]
        similarities = np.array(atom1) - np.array(atom2)
        if similarities[0] == 0:
            intra_mon += 1
        elif similarities[1] == 0:
            intra_layer += 1
        elif similarities[2] == 0:
            intra_pore += 1
        else:
            inter_pore += 1

    # Output

    print('')
    print('Total Number of crosslinks: %s' % nxlinks)
    print(' ______________________________________')
    print('|____________________| Number |   %    |')
    print('|Intra-Monomer xlinks|{:^8}|{:^8.2f}|'.format(intra_mon, 100*intra_mon/nxlinks))
    print('|Intra-Layer xlinks  |{:^8}|{:^8.2f}|'.format(intra_layer, 100*intra_layer/nxlinks))
    print('| Intra-Pore xlinks  |{:^8}|{:^8.2f}|'.format(intra_pore, 100*intra_pore/nxlinks))
    print('| Inter-Pore xlinks  |{:^8}|{:^8.2f}|'.format(inter_pore, 100*inter_pore/nxlinks))
    print(' --------------------------------------')