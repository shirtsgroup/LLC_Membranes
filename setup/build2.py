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
    parser.add_argument('-r', '--pore_radius', default=6, type=float, help='Initial Pore Radius (Angstroms)')
    parser.add_argument('-p', '--p2p', default=45, type=float, help='Initial Pore to Pore Distance')
    parser.add_argument('-n', '--nopores', default=4, type=int, help='Number of Pores')
    parser.add_argument('-d', '--dbwl', default=3.7, type=float, help='Distance between layers')
    parser.add_argument('-t', '--tilt', default=0, type=float, help='Monomer tilt angle')
    parser.add_argument('-O', '--offset', help="Specify this flag to build the system in an offset configuration",
                        action="store_true")
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
    parser.add_argument('-L', '--correlation_length', default=10, type=float, help='Desired correlation length')
    parser.add_argument('-Lvar', default=0.1, type=float, help='variance in z position of monomer heads')

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

    def __init__(self, name, npores, p2p, pore_alpha, ncolumns, monomers_per_column):
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

        self.pore_centers = np.zeros([npores, 2])

        pore_alpha_radians = pore_alpha * (np.pi / 180)
        self.pore_centers[1, :] = [p2p*np.cos(pore_alpha_radians), p2p*np.sin(pore_alpha_radians)]
        self.pore_centers[2, :] = [p2p*np.cos(pore_alpha_radians) + p2p, p2p*np.sin(pore_alpha_radians)]
        self.pore_centers[3, :] = [p2p, 0]

    def build_column(self):


if __name__ == "__main__":

    args = initialize().parse_args()

    system = Assembly(args.build_monomer, args.nopores, args.p2p, args.angles[2], args.ncolumns,
                      args.monomers_per_column)

    for i in range(args.nopores):
        for j in range(args.ncolumns):
            z = np.linspace(0, args.dbwl*args.monomers_per_column - args.dbwl, args.monomers_per_column)
            print(z)
            exit()