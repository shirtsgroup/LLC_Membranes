#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from llclib import physical


def initialize():

    parser = argparse.ArgumentParser(description='Measure the density inside the hydrophobic region of an LLC system')

    parser.add_argument('-t', '--traj', default='wiggle.xtc', type=str, help = 'Trajectory file (.xtc or .trr). The'
                        'trajectory should be preprocessed beforehand to make molecules whole or else the distance '
                        'calculation get thrown off. Use "gmx trjconv -pbc whole" to do this')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help = 'Name of coordinate file')
    parser.add_argument('-c', '--components', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'],
                        help='component used to track pore positions and define pore radius')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)  # load gromacs trajectory
    keep = [a.index for a in t.topology.atoms if a.name in args.components]   # get the index of all atoms in components
    pore_components = t.atom_slice(keep)  # create a new trajectory object only describing those atoms
    pos = pore_components.xyz

    pcenters = physical.avg_pore_loc(4, pos)  # find the pore centers

    r, r_std = physical.limits(pos, pcenters)

    print 'Average Pore Size: %1.2f nm' % np.mean(r)