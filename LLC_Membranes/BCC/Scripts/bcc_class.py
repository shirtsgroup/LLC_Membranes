#!/usr/bin/env python

import os
# from LLC_Membranes.analysis import Atom_props
import numpy as np

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


class LC(object):
    """A Liquid Crystal monomer has the following attributes which are relevant to building and crosslinking:

    Attributes:
        name: A string representing the monomer's name.
        natoms: An integer accounting for the number of atoms in a single monomer.
        build_mon: Monomer used to build the unit cell
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
        with open('%s/../top/structures/%s' % (location, name)) as f:
            for line in f:
                a.append(line)

        P = []
        L1 = []  # line group 1
        L2 = []  # line group 2
        R = []
        R_i = []  # reference indices
        self.no_ions = 0
        self.ions = []
        names = []
        resid = []
        MW = 0  # molecular weight
        MW = 716

        if name.endswith('.gro'):

            # Lists of Plane atoms, line atoms, and dummy atoms
            self.natoms = int(a[1])
            xyz = np.zeros([self.natoms, 3])
            residues = [str.strip(a[2][5:10])]
            nres = []
            res_count = 0
            self.ref_names = []
            for i in range(2, len(a) - 1):
                xyz[i - 2, :] = [float(a[i][20:28]), float(a[i][28:36]), float(a[i][36:44])]
                res = str.strip(a[i][5:10])
                name = str.strip(a[i][10:15])
                #MW += Atom_props.mass[name]
                names.append(name)
                resid.append(res)
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
                        P.append([float(a[i][20:28]), float(a[i][28:36]), float(a[i][36:44])])
                    if 'L1' in annotations:
                        L1.append([float(a[i][20:28]), float(a[i][28:36]), float(a[i][36:44])])
                    if 'L2' in annotations:
                        L2.append([float(a[i][20:28]), float(a[i][28:36]), float(a[i][36:44])])
                    if 'R' in annotations:
                        R.append([float(a[i][20:28]), float(a[i][28:36]), float(a[i][36:44])])
                        R_i.append(i - 2)
                        self.ref_names.append(str.strip(a[i][10:15]))
                    if 'I' in annotations:
                        self.no_ions += 1
                        ion = str.strip(a[i][10:15])
                        #self.valence = Atom_props.charge[ion]
                        if ion not in self.ions:
                            self.ions.append(ion)

            nres.append(res_count)

        elif name.endswith('.pdb'):

            # Lists of Plane atoms, line atoms, and dummy atoms

            start = 0
            while a[start].count('ATOM') == 0:
                start += 1
            end = start
            while a[end].count('ATOM') != 0:
                end += 1

            self.natoms = end - start
            xyz = np.zeros([self.natoms, 3])

            residues = [str.strip(a[start][17:22])]
            nres = []
            res_count = 0
            for i in range(start, end):
                line = a[i].split()
                xyz[i - start, :] = [float(line[5]), float(line[6]), float(line[7])]
                res = str.strip(a[i][17:22])
                name = str.strip(a[i][11:17])
                #MW += Atom_props.mass[name]
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
                    if 'L1' in annotations:
                        L1.append([float(line[5]), float(line[6]), float(line[7])])
                    if 'L2' in annotations:
                        L2.append([float(line[5]), float(line[6]), float(line[7])])
                    if 'R' in annotations:
                        R.append([float(line[5]), float(line[6]), float(line[7])])
                    if 'I' in annotations:
                        self.no_ions += 1
                        ion = line[2]
                        #self.valence = Atom_props.charge[ion]
                        if ion not in self.ions:
                            self.ions.append(ion)

            nres.append(res_count)

        # make everything into numpy arrays
        Plane = np.zeros([3, 3])
        for i in range(3):
            Plane[i, :] = P[i]
	
        # direction monomer is pointing
        L = np.zeros([2, 3])
        L[0, :] = np.array([sum([l[0] for l in L1]), sum([l[1] for l in L1]), sum([l[2] for l in L1])]) / len(L1)
        L[1, :] = np.array([sum([l[0] for l in L2]), sum([l[1] for l in L2]), sum([l[2] for l in L2])]) / len(L2)
        V = L[1, :] - L[0, :]  # direction monomer is pointing based on points defined by average xyz of L1 and L2

        # reference position
        ref = np.array([sum([r[0] for r in R]), sum([r[1] for r in R]), sum([r[2] for r in R])]) / len(R)

        self.xyz = xyz
        self.planeatoms = P
        self.linevector = V
        self.reference = ref
        self.residues = residues
        self.nresidues = nres
        self.ref_index = R_i
        self.resid = resid
        self.names = names
        self.linepts = L
        self.MW = MW
