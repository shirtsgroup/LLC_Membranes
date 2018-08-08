#! /usr/bin/env python3

import numpy as np
import argparse
from LLC_Membranes.llclib import file_rw, transform
from LLC_Membranes.setup.lc_class import LC
import os
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Build HII LLC unit cell')

    parser.add_argument('-b', '--build_monomer', default='NAcarb11V.gro', type=str, help='Name of class of monomer using to build with')
    parser.add_argument('-o', '--out', default='initial.gro', help='name of output file')
    parser.add_argument('-m', '--monomers_per_column', default=20, type=int, help='Number of monomers per column')
    parser.add_argument('-c', '--ncolumns', default=5, type=int, help='Number of columns per pore')
    parser.add_argument('-r', '--pore_radius', default=.6, type=float, help='Initial Pore Radius (nm)')
    parser.add_argument('-p', '--p2p', default=4.5, type=float, help='Initial Pore to Pore Distance (nm)')
    parser.add_argument('-n', '--nopores', default=4, type=int, help='Number of Pores')
    parser.add_argument('-d', '--dbwl', default=.37, type=float, help='Distance between layers (nm)')
    parser.add_argument('-t', '--tilt', default=0, type=float, help='Monomer tilt angle')
    parser.add_argument('-O', '--offset', help="Specify this flag to build the system in an offset configuration",
                        action="store_true")
    parser.add_argument('-L', '--correlation_length', type=float, help='Desired correlation length')
    parser.add_argument('-Lvar', default=0.1, type=float, help='variance in z position of monomer heads')
    parser.add_argument('-pd', '--parallel_displaced', default=0, type=float, help='Distance between center of mass '
                        'of vertically stacked monomers on the xy plane. A distance of 0 is the same as sandwiched.')
    parser.add_argument('--rot', default=45, type=float, help="Rotate pores by this amount (degrees)")
    parser.add_argument('--offset_angle', default=0, type=float)
    parser.add_argument('-box', '--box_lengths', nargs='+', type=float, help='box vector lengths [x y z]')
    parser.add_argument('-angles', '--angles', nargs='+', default=[90, 90, 60], type=float, help='angles between box'
                                                                                                  'box vectors')
    parser.add_argument('-rd', '--radial_displacement', type=float, help='Shift pore center by this value for every'
                                                                         'other layer')
    parser.add_argument('-rdtheta', default=0, type=float, help='When radially displaced, the angle with respect'
                        'to the x-axis defining the direction in which to shift the pore center')
    parser.add_argument('-ad', '--angularly_displaced', action="store_true", help='rotate phenyl groups'
                        'with respect to those in adjacent layers')
    parser.add_argument('-columns', action="store_true", help='build column-wise')


    return parser


def z_correlation(z, L, v=0.1):
    """
    Calculate where to place monomers on the z-axis so that a given correlation length is obtained
    :param z: mean z-positions where monomers will be placed with gaussian probability np.array([n_layers])
    :param L: desired correlation length [float]
    :param v: variance in z position of monomer head groups
    :return: locations [np.array[nlayers])
    """

    n = z.shape[0]
    cov = np.zeros([n, n])  # initialize covariance matrix

    decay = v*np.exp(-z / L)  # decay of covariance
    # decay[1:] += np.exp(-z[::-1][:-1]/L) # for periodicity (?)

    for i in range(z.shape[0]):
        cov[i, i:] = decay[:(n - i)]
        cov[i:, i] = decay[:(n - i)]

    locations = np.random.multivariate_normal(z, cov)

    return locations


class Assembly(LC):

    def __init__(self, name, npores, p2p, pore_alpha, pore_radius, tilt=0):

        """Initialize geometry of columnar pore structure

        Keyword arguments:
            name -- name of monomer with which the system will be built
            npores -- number of pores in the system
            p2p -- absolute pore-to-pore distance
            pore_alpha -- angle between x and y box vector. For example if pore_alpha = 120 or 60, you'll get
            hexagonally packed pores
            ncolumns -- number of columns surrounding each pore
            monomers_per_column -- number of monomers in each column
        """

        super().__init__(name)

        self.tilt = tilt
        self.xyz = np.zeros([0, 3])
        self.names = []
        self.all_residues = []
        self.pore_radius = pore_radius

        # currently only implemented for 4 pores
        self.pore_centers = np.zeros([npores, 2])

        pore_alpha_radians = pore_alpha * (np.pi / 180)

        self.pore_centers[1, :] = [p2p*np.cos(pore_alpha_radians), p2p*np.sin(pore_alpha_radians)]
        self.pore_centers[2, :] = [p2p*np.cos(pore_alpha_radians) + p2p, p2p*np.sin(pore_alpha_radians)]
        self.pore_centers[3, :] = [p2p, 0]

        # center pores in box
        self.pore_centers[:, 0] += (p2p / 2) * (1 + np.cos(pore_alpha_radians))
        self.pore_centers[:, 1] += (p2p / 2) * np.sin(pore_alpha_radians)

    def align_plane(self):
        """ Align the atoms defined by the plane_indices attribute of LC with the xy plane """

        plane_atoms = np.zeros([3, 3])
        for i in range(plane_atoms.shape[0]):
            plane_atoms[i, :] = self.LC_positions[self.plane_indices[i], :]

        R = transform.rotateplane(plane_atoms, angle=self.tilt)  # generate rotation matrix

        b = np.ones([1])
        for i in range(self.LC_positions.shape[0]):
            coord = np.concatenate((self.LC_positions[i, :], b))
            x = np.dot(R, coord)
            self.LC_positions[i, :] = x[:3]

    def translate_to_origin(self):
        """ Translate molecule to the origin using the ref_atom_index attribute of LC """

        self.LC_positions = transform.translate(self.LC_positions, self.LC_positions[self.ref_atom_index, :], [0, 0, 0])

    def align_with_x(self):
        """ Align vector defined by lineatoms in LC object with x axis """

        v = np.array([self.LC_positions[self.lineatoms[0], :2] - self.LC_positions[self.lineatoms[1], :2]])
        angle = np.arctan(v[0, 1] / v[0, 0])
        self.LC_positions = transform.rotate_coords_z(self.LC_positions, - angle * 180 / np.pi)

    def build_column(self, pore, z, theta, correlation=True, var=0, correlation_length=0, pd=0, random_shift=True):
        """ Place a column at angle theta on xy plane with respect to a pore center

        Keyword Arguments:
            pore -- pore number (0 : npores - 1)
            z -- mean z-positions of monomers in column
            theta -- angle, with respect to pore center where column should be placed
            correlation -- adjust z positions so there is a correlation length
            var -- variance in multivariate normal distribution used to make correlated points
            correlation_length -- length for which correlation between stacked monomers to persist
            pd -- Specify a nonzero value to make a parallel displaced. The distance here (in nm) will be the distance
            between the center of masses of two vertically stacked monomers.
        """

        if correlation:
            z = z_correlation(z, correlation_length, v=var)
            if random_shift:
                dbwl = z[1] - z[0]  # distance between stacked monomers
                z += np.random.uniform(0, dbwl)
        elif random_shift:
            dbwl = z[1] - z[0]  # distance between stacked monomers
            z += np.random.uniform(0, dbwl)

        displaced_theta = (180 / np.pi) * (2 * np.arcsin(pd / (2 * self.pore_radius)))

        natoms = self.LC_positions.shape[0]  # number of atoms including ions
        pos = np.copy(self.LC_positions)
        pos[:, 0] += self.pore_radius
        displaced_pos = np.copy(pos)

        pos = transform.rotate_coords_z(pos, theta)
        displaced_pos = transform.rotate_coords_z(displaced_pos, theta + displaced_theta)

        column = np.zeros([z.size * natoms, 3])

        before = np.array([0, 0, 0])
        for l in range(z.size):
            if l % 2 == 0:
                column[l * natoms:(l + 1) * natoms, :] = transform.translate(pos, before, np.array([0, 0, z[l]]))
            else:
                column[l * natoms:(l + 1) * natoms, :] = transform.translate(displaced_pos, before, np.array([0, 0, z[l]]))
            self.names += self.LC_names
            self.all_residues += self.LC_residues

        column = transform.translate(column, before, [self.pore_centers[pore, 0], self.pore_centers[pore, 1], 0])

        self.xyz = np.concatenate((self.xyz, column))

    def reorder(self):
        """ reorder coordinate, residues and atom names so that residues are separated """

        residue_indices = []
        for r in self.residues:
            ndx = [i for i, a in enumerate(self.all_residues) if a == r]
            residue_indices.append(ndx)

        ordered = []
        for i in range(len(residue_indices)):
            ordered += residue_indices[i]

        self.xyz = self.xyz[ordered, :]
        self.all_residues = [self.all_residues[i] for i in ordered]
        self.names = [self.names[i] for i in ordered]

    def write_gro(self, out, ucell):

        file_rw.write_gro_pos(self.xyz, out, ids=self.names, res=self.all_residues, ucell=ucell)


if __name__ == "__main__":

    args = initialize().parse_args()

    correlation = False
    if args.correlation_length is not None:
        correlation = True

    system = Assembly(args.build_monomer, args.nopores, args.p2p, args.angles[2], args.pore_radius)

    system.align_plane()  # align monomer head group with xy plane
    system.translate_to_origin()  # move monomer to origin for rotation
    system.align_with_x()  # align vector from benzene ring to carboxylate with x axis

    wedge_theta = 360 / args.ncolumns  # rotation between laterally adjacent monomers (angle defining slice)
    for i in range(args.nopores):
        start_theta = np.random.uniform(0, 360)
        thetas = [start_theta + x*wedge_theta for x in range(args.ncolumns)]
        for j in range(args.ncolumns):
            z = np.linspace(0, args.dbwl*args.monomers_per_column - args.dbwl, args.monomers_per_column)
            system.build_column(i, z, thetas[j], correlation=correlation, var=args.Lvar,
                                correlation_length=args.correlation_length, pd=args.parallel_displaced, random_shift=True)

    system.reorder()

    if args.box_lengths:
        a, b, c = args.box_lengths
    else:
        a, b, c = [2*args.p2p, 2*args.p2p, args.dbwl*args.monomers_per_column]  # logical choices for symmetry

    alpha, beta, gamma = [angle * (np.pi / 180) for angle in args.angles]

    V = a * b * c * np.sqrt(
        1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.sin(
            gamma))  # volume of unitcell

    A = np.array([a, 0, 0])  # vector in x direction
    B = np.array([b*np.cos(gamma), b*np.sin(gamma), 0])  # vector in y direction
    C = np.array([c*np.cos(beta), c*((np.cos(alpha) - np.cos(gamma)*np.cos(beta))/np.sin(gamma)), V / (a*b*np.sin(gamma))])  # vector in z direction

    box = np.vstack((A, B, C))

    system.write_gro(args.out, box)
