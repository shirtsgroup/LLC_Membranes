#!/usr/bin/env python

import os
from LLC_Membranes.analysis import Atom_props
import mdtraj as md
import numpy as np

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


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
        with open('%s/../top/topologies/%s' % (location, name)) as f:
            for line in f:
                a.append(line)

        t = md.load("%s/../top/topologies/%s" % (location, name))
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
        self.no_ions = 0
        self.ions = []
        self.MW = 0

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

# test = LC('NAcarb11V.pdb')
# print test.natoms
# print test.planeatoms
# print test.lineatoms
# print test.ref_atom_index
# print test.residues
# print test.nresidues
# print test.c1_atoms
# print test.c2_atoms
# print test.valence
