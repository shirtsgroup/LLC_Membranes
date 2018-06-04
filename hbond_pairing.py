#! /usr/bin/env python

"""
Calculate H-bonds in system.

Some definitions:
    donor (D): atom covalently bonded to hydrogen
    acceptor (A) : atom which 'accepts' hydrogen bond

A diagram of an hbond:

    D--H - - A

Criterion:
    Distance between D and A below some distance
    Angle between DHA less than some cut-off

Exact values for criterion are left up to the user.
"""


import argparse
import mdtraj as md
import numpy as np
import os
import tqdm
import matplotlib.pyplot as plt
import hbonds

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='wiggle.trr', type=str, help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate file')
    parser.add_argument('-p', '--top', default='topol.top', type=str, help='Gromacs topology file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='End frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Start frame')
    parser.add_argument('-x', '--exclude_water', action='store_true', help='Exclude water while searching for hbonds')
    parser.add_argument('-r', '--residues', nargs='+', default=['HII'], help='Residues to include in h-bond search')
    parser.add_argument('-a', '--atoms', action='append', nargs='+', help='Atoms for to include for each '
                        'residue. Each list of atoms must be passed with a separate -a flag for each residue in'
                        'args.residues')
    parser.add_argument('-d', '--distance', default=.3, help='Maximum distance between acceptor and donor atoms')
    parser.add_argument('-angle', '--angle_cut', default=20, help='Maximum DHA angle to be considered an H-bond')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    # workaround for argparse. If default value is set, it is always included in the list with action='append'
    if not args.atoms:
        args.atoms = ['O3', 'O4']  # a default value

    sys = hbonds.System(args.traj, args.gro, args.top, begin=args.begin, end=args.end)

    for i, r in enumerate(args.residues):
        sys.set_eligible(r, args.atoms[i])

    sys.identify_hbonds(args.distance, args.angle_cut)
    sys.plot_hbonds()