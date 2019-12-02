#!/usr/bin/env python

"""
    The purpose of this script is to edit the topology file of a system containing molecules which have a benzene ring
    in order to create your choice of two things:
    (1) An artificial dipole which will act as electron clouds participating in a pi bond. The dipole is
        created by centering two virtual sites above and below the plane of the benzene ring and assigning them
        appropriate charges values.
    (2) Add position restraints with a given force constant to chosen atoms w.r.t. to a specified axis or axes
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import argparse
import numpy as np
import warnings
import os
from LLC_Membranes.llclib import file_rw, topology
import mdtraj as md

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Location of this script


def initialize():

    parser = argparse.ArgumentParser(description='Duplicate points periodically in the x-y directions')

    parser.add_argument('-g', '--gro', default='initial.gro', type=str, help='Coordinate file')
    parser.add_argument('-o', '--out', default='dipole.itp', type=str, help='Name of output topology file')
    parser.add_argument('-f', '--f_const', nargs='+', default=[1000, 1000, 1000], type=float, help='Force constant')
    parser.add_argument('-a', '--atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], type=str, help='Name of carbons in ring')
    parser.add_argument('-d', '--distance', default=0.1, help='Distance to offset dipole from ring (Angstroms)')
    parser.add_argument('-m', '--monomer', default='NAcarb11V', help='Which monomer topology is being used')
    parser.add_argument('-A', '--axis', default='xy', help='Axis to restrain along with position restraints')
    ######## These parameters are not implemented yet but left here since they will have some use in the future ########
    parser.add_argument('--xlink', help='Specify this flag if the system is being crosslinked while restrained',
                        action="store_true")
    parser.add_argument('--append', help='Specify this to prevent the topology from being re-written, and instead, add '
                                         'restraints to the topology for other atoms', action="store_true")
    parser.add_argument('-i', '--input', type=str, default='dipole.itp', help='Name of topology file to be edited if '
                                                                              'you use the option --append')
    ####################################################################################################################
    parser.add_argument('-dr', '--dihedral_restraints', action='append', nargs='+',
                        help='Specify atom names of dihedral to be restrained, followed by angle at which to restrain'
                             'them, the deviation from that angle allowed and the force constant to apply. For example:'
                             '"restrain.py -dr C1 C C6 O4 90 0 1000" means to keep the angle between the planes formed'
                             'by C1-C-C6 and C-C6-O4 90 degrees apart with 0 degrees of leeway and a force constant of '
                             '1000')
    parser.add_argument('-com', '--center_of_mass', action="store_true", help="Add position restraints at the center of"
                                                                              "mass of args.atoms")
    parser.add_argument('-v', '--virtual_site_parameters', nargs='+', default=['3fd', 'C', 'C2', 'C4', '.5', '.14'],
                        help='A list in the following order : virtual site construction type, atoms to use to build '
                             'virtual site, required length parameters (i.e. a, b, c or d) in the order specified in'
                             'GROMACS documentation')
    parser.add_argument('-b', '--bond_restraints', default=[.37, 1000], nargs='+', help='Bond restraint pararmeters. A'
                        'list where the first entry is the equilibrium distance (nm) and the second entry is the force'
                        'constant for a harmonic potential (kJ/mol/nm^2)')

    args = parser.parse_args()

    return args


warnings.filterwarnings("error")  # This makes it so numpy warnings are treated as real errors


def virtual_sites(all_coords, monomers, valence, a, b, c, funct):

    n_atoms = np.shape(all_coords)[1] - (1 / valence) * monomers
    atoms_per_molecule = n_atoms / monomers  # subtract valence to exclude counterion

    n_vsites = monomers * 2  # number of virtual sites need to create a dipole (2 per ring)
    vsites = np.zeros([8, n_vsites])
    vsites[4, :] = funct
    vsites[5, :] = a  # all a and b values are the same
    vsites[6, :] = b

    for i in range(monomers):
        for j in range(2):
            vsites[0, i * 2 + j] = n_atoms + i * 2 + j + 1  # new atoms were placed at the end of the .gro file
            vsites[1:4, i * 2 + j] = [i * atoms_per_molecule + 1, i * atoms_per_molecule + 3,
                              i * atoms_per_molecule + 5]
            vsites[7, i * 2 + j] = (-1)**j * c

    return vsites


def exclusions(coord_file, monomers, valence, toplines, atoms, n_atoms, vsites):
    """
    :param coord_file: the original .gro file stored in a list
    :param monomers: number of monomers
    :param valence: charge on ions
    :param atoms: the names of the atoms which should be excluded
    :param n_atoms: the number of atoms total
    :return: a list of exclusions
    """
    n_atoms = np.shape(all_coords)[1] - (1 / valence) * monomers
    atoms_per_molecule = n_atoms / monomers  # subtract valence to exclude counterion (includes new PI atoms)

    n_excluded = len(atoms) + 1  # number of atoms excluded. +1 because it will be excluded from it complementary vsite

    n_exclusions = monomers * 2  # number of exclusions to be specified (one for each vsite)
    exclusions = np.zeros([n_excluded + 1, n_exclusions])  # the first entry is the virtual site itself

    for i in range(np.shape(vsites)[1]):
        exclusions[0, i] = vsites[0, i]
        # add exclusion from the complementary vsite. This is pretty specific to the format and will likely need to be
        # re-written eventually
        if i % 2 == 0:
            exclusions[1, i] = vsites[0, i + 1]
        elif i % 2 == 1:
            exclusions[1, i] = vsites[0, i - 1]

    x = 0
    for i in range(monomers):
        a = 2
        for j in range(atoms_per_molecule):
            line = i*atoms_per_molecule + j + toplines
            if str.strip(coord_file[line][10:15]) in atoms:
                exclusions[a, x] = int(coord_file[line][15:20])
                exclusions[a, x + 1] = int(coord_file[line][15:20])
                a += 1
        x += 2

    return exclusions


def dihedral_restraints(file, atoms):
    """
    This function needs to be moved into the class
    """

    ndihedrals = len(atoms)

    all_restraints = np.zeros([0, 8])
    for n in range(ndihedrals):

        atom_numbers = []
        d = np.zeros([4])
        count = 0
        for line in file:
            atom = str.strip(line[10:15])  # name of atom at that line
            if atom in atoms[n]:
                d[atoms[n].index(atom)] = int(line[15:20])  # atom number placed in d in order that dihedral was passed
                if np.count_nonzero(d) == 4:  # and len(atom_numbers) % 8 == 0:
                    if count % 2 == 1:
                        for i in range(4):
                            atom_numbers.append(d[i])
                    d = np.zeros([4])
                    count += 1

        restraints = np.zeros([len(atom_numbers)//4, 8])
        for i in range(len(atom_numbers)//4):
            d = 4 * i
            restraints[i, :] = [atom_numbers[d], atom_numbers[d + 1], atom_numbers[d + 2], atom_numbers[d + 3], 1, atoms[n][4],
                                atoms[n][5], atoms[n][6]]

        all_restraints = np.concatenate((all_restraints, restraints))

    return all_restraints


class RestrainedTopology(object):

    def __init__(self, gro, res, atoms, name='restrained', com=False, xlink=False,
                 vparams=None):
        """ Write topology to restrain one or more residues with position restraints in GROMACS

        :param gro: coordinate file where restraints will be placed
        :param res: name of residue where position restraints are being added
        :param atoms: name of atoms to be restrained in res
        :param name: name of output topology file
        :param com: restrain center of mass of atoms instead of individual atoms
        :param xlink : whether or not the system is in the process of being crosslinked
        :param vparams: A list in the following order : virtual site construction type, atoms to use to build virtual
               site, required length parameters (i.e. a, b, c or d) in the order specified in GROMACS documentation
        """

        topology.fix_resnumbers(gro)
        t = md.load(gro)

        if type(res) is str:
            res = [res]
            atoms = [atoms]

        self.all_coords = t.xyz[0, :, :]  # all coordinates for system
        self.atom_numbers = []
        for i in range(len(atoms)):
            self.atom_numbers.append([a.index + 1 for a in t.topology.atoms if a.name in atoms[i] and a.residue.name ==
                                      res[i]])

        self.atoms = t.n_atoms  # number of atoms in full system
        self.nmon = [len(self.atom_numbers[i]) // len(atoms[i]) for i in range(len(atoms))]  # number of monomer residues
        self.name = name  # name of output files (.itp, .gro if you are using centers of masses)
        self.residue = res  # name of residue(s) to which position restraints are being applied
        self.LC = [topology.LC('%s' % r) for r in self.residue]  # everything we can know about the residue
        self.com = com

        # self.keep = np.array([a.index for a in t.topology.atoms if a.name in atoms])  # atoms to keep
        # self.atom_numbers = self.keep + 1  # numbers (not indices) of atoms to keep
        # self.coords = self.all_coords[self.keep, :]  # coordinates of atoms in keep

        if self.com:  # add center of mass virtual site

            # These things are only needed for center of mass virtual site construction
            self.ids = [a.name for a in t.topology.atoms]  # names of all atoms in system
            self.res = [a.residue.name for a in t.topology.atoms]  # residue names of all atoms in system
            self.vparams = vparams
            self.vatoms_numbers = [a.index + 1 for a in t.topology.atoms if a.name in vparams]
            self.box = t.unitcell_vectors[0, :, :]  # box vectors in mdtraj formate
            self.box_gromacs = [self.box[0, 0], self.box[1, 1], self.box[2, 2], self.box[0, 1], self.box[2, 0],
                                self.box[1, 0], self.box[0, 2], self.box[1, 2],
                                self.box[2, 0]]  # gromacs format box vects

            with open('%s/../top/Monomer_Tops/%s.itp' % (location, self.residue), 'r') as f:
                residue_top = []
                for line in f:
                    residue_top.append(line)

            atoms_index = 0
            while residue_top[atoms_index].count('[ atoms ]') == 0:
                atoms_index += 1

            while residue_top[atoms_index] != '\n':
                atoms_index += 1

            residue_top.insert(atoms_index, '{:>6d}{:>5s}{:>6d}{:>6s}{:>6s}{:>5d}{:>13.6f}{:>13.6f}\n'.format(
                self.LC.natoms + 1, 'hc_d', 1, self.LC.residues[0], 'HD', self.LC.natoms + 1, 0, 0))

            if self.vparams[0] == '3fd':
                # 'a' = 0.5 and 'b' = 0.14 (aromatic carbon bond length (nm)) puts a vsite in the middle of benzene
                # if the constructor atoms are 3 non-adjacent carbons from the ring
                residue_top.append('[ virtual_sites3 ]\n')
                residue_top.append('{:<6d}{:<6d}{:<6d}{:<6d}{:<6d}{:<8.4f}{:<8.4f}\n'.format(self.LC.natoms + 1,
                                    self.vatoms_numbers[0], self.vatoms_numbers[1], self.vatoms_numbers[2], 2,
                                    float(self.vparams[-2]), float(self.vparams[-1])))
            else:
                print('Your choice of virtual site has not yet been implemented')
                exit()

            # groups = np.reshape(self.coords, (len(self.keep) // len(atoms), len(atoms), 3))
            # centers_of_mass = np.mean(groups, axis=1)
            # self.coords = centers_of_mass  # redefine coordinates as centers of mass

            file_rw.write_assembly(residue_top, '%s.itp' % self.name, self.nmon, xlink=xlink)

            # now the dummies need to be added to the .gro file. They are placed at the end of the residue section
            # This loop works for a single virtual site per monomer. It will need to be modified if multiple sites
            # are to be constructed.
            insert_ndx = self.LC.natoms
            self.atom_numbers = []  # redefine this since everything is renumbered
            for i in range(self.nmon):
                ndx = (i + 1)*insert_ndx + i
                self.ids.insert(ndx, 'HD')
                self.res.insert(ndx, 'HII')  # should make this more general
                self.all_coords = np.insert(self.all_coords, ndx, np.array([0, 0, 0]), axis=0)
                self.atom_numbers.append(ndx + 1)

            file_rw.write_gro_pos(self.all_coords, '%s.gro' % self.name, ids=self.ids, res=self.res, box=self.box_gromacs)
        else:

            file_rw.write_assembly(res, '%s.itp' % self.name, self.nmon, xlink=xlink)

        with open('%s.itp' % self.name, 'r') as f:
            self.topology = []
            for line in f:
                self.topology.append(line)

    def add_position_restraints(self, axis, f_const):
        """
        Restrain the selected atoms in desired directions
        :param axis: which direction to restrain (xyz, xy, z, xz .. etc.)
        :param f_const: force constant in each direction. Order of force constants matches that of axis argument
        :return: an array of position restraints formatted for easy writing into the topology (.itp)
        """

        self.topology.append("\n[ position_restraints ]\n")

        fc = np.zeros([3])
        for i, a in enumerate(axis):
            if a == 'x':
                fc[0] = f_const[i]
            if a == 'y':
                fc[1] = f_const[i]
            if a == 'z':
                fc[2] = f_const[i]

        atom_numbers = []
        for i in self.atom_numbers:
            atom_numbers += i

        restraints = np.zeros([5, len(atom_numbers)])  # organize them into a list which can be translated to a topology
        for i in range(len(atom_numbers)):
            restraints[:, i] = [atom_numbers[i], 1, fc[0], fc[1], fc[2]]  # See: http://www.gromacs.org/Documentation/How-tos/Position_Restraints
            self.topology.append('{:6d}{:6d}{:1s}{:9f}{:1s}{:9f}{:1s}{:9f}\n'.format(int(restraints[0, i]),
                                int(restraints[1, i]),'', restraints[2, i], '', restraints[3, i], '', restraints[4, i]))

    def add_distance_restraint_columns(self, b0, kb, layers=20, pores=4):
        """
        Add distance constraints to centers of mass of monomer head groups. This is a function specialized for an
        HII system built with build.py (without the flag -columns).
        :param b0 : equilibrium distance
        :param kb : force constant for harmonic potential
        :param layers : layers per pore
        :param pores : number of pore columns
        """

        self.topology.append("\n[ bonds ]\n")
        mpl = int(self.nmon / 4 / layers)  # monomers per layer

        for p in range(pores):
            for l in range(layers):
                for m in range(mpl):
                    if l == (layers - 1):  # handles periodicity
                        self.topology.append('{:<6d}{:<6d}{:<6d}{:<6.1f}{:<6.1f}\n'.format(self.atom_numbers[layers*mpl*p
                                             + l*mpl + m], self.atom_numbers[layers*mpl*p + m], 6, b0, kb))
                    else:
                        self.topology.append('{:<6d}{:<6d}{:<6d}{:<6.1f}{:<6.1f}\n'.format(self.atom_numbers[layers*mpl*p
                                             + l*mpl + m], self.atom_numbers[layers*mpl*p + (l + 1)*mpl + m], 6, b0, kb))

        # Tether together columns (this is temporary so variables are hard coded)
        b0 = 2*0.6*np.sin(36*np.pi/180)
        kb = 1000
        for p in range(pores):
            for l in range(layers):
                for m in range(mpl):
                    if m == (mpl - 1):  # handles periodicity
                        self.topology.append('{:<6d}{:<6d}{:<6d}{:<6.3f}{:<6.1f}\n'.format(self.atom_numbers[layers*mpl*p
                                             + l*mpl + m], self.atom_numbers[layers*mpl*p + l*mpl], 6, b0, kb))
                    else:
                        self.topology.append('{:<6d}{:<6d}{:<6d}{:<6.3f}{:<6.1f}\n'.format(self.atom_numbers[layers*mpl*p
                                             + l*mpl + m], self.atom_numbers[layers*mpl*p + l*mpl + m + 1], 6, b0, kb))

    def write_topology(self):
        with open('%s.itp' % self.name, 'w') as f:
            for line in self.topology:
                f.write(line)


if __name__ == "__main__":

    args = initialize()

    top = RestrainedTopology(args.gro, args.monomer, args.atoms, com=args.center_of_mass,
                             vparams=args.virtual_site_parameters)
    # top.add_distance_restraint_columns(float(args.bond_restraints[0]), float(args.bond_restraints[1]))
    top.add_position_restraints(args.axis, args.f_const)
    top.write_topology()
