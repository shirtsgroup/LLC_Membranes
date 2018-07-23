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
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of Layers')
    parser.add_argument('-m', '--monomers', default=5, type=int, help='Monomers per layer')
    parser.add_argument('-r', '--radius', default=6, type=float, help='Initial Pore Radius (Angstroms)')
    parser.add_argument('-p', '--p2p', default=45, type=float, help='Initial Pore to Pore Distance')
    parser.add_argument('-n', '--nopores', default=4, type=int, help='Number of Pores')
    parser.add_argument('-d', '--dbwl', default=3.7, type=float, help='Distance between layers')
    parser.add_argument('-s', '--layer_distribution', default='uniform', help='The distribution of monomers per layer')
    parser.add_argument('-a', '--alt_1', default=6, type=int, help='Monomers per layer for the first type of alternating layer')
    parser.add_argument('-A', '--alt_2', default=8, type=int, help='Monomers per layer for the second type of alternating layer')
    parser.add_argument('-t', '--tilt', default=0, type=float, help='Monomer tilt angle')
    parser.add_argument('-H', '--helix', help="Specify this flag if you want to build in a helical configuration",
                        action="store_true")
    parser.add_argument('-O', '--offset', help="Specify this flag to build the system in an offset configuration",
                        action="store_true")
    parser.add_argument('--rot', default=45, type=float, help="Rotate pores by this amount (degrees)")
    parser.add_argument('--offset_angle', default=0, type=float)
    parser.add_argument('-box', '--box_lengths', nargs='+', type=float, help='box vector lengths '
                                                                                                '[x y z] ')
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


class Assembly(LC):

    def init(self, name):

        super().__init__(name)


if __name__ == "__main__":

    args = initialize().parse_args()

    system = Assembly(args.build_monomer)