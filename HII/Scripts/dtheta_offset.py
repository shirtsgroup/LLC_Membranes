#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from llclib import physical
import math


def initialize():

    parser = argparse.ArgumentParser(description='Measure Membrane Thickness')

    parser.add_argument('-g', '--gro', type=str, default='wiggle.gro', help='Path to input file')
    parser.add_argument('-t', '--traj', type=str, default='wiggle.trr', help='Trajectory file (.xtc or .trr)')
    parser.add_argument('-r', '--ref', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='reference atoms for '
                                                                'calculating pore centers and radii')
    parser.add_argument('-a', '--layer_offset', type=float, help='angle of offset between layers')
    parser.add_argument('-d', '--dbwl', default=.37, type=float, help='distance between layers (nm)')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)
    keep = [a.index for a in t.topology.atoms if a.name in args.ref]
    pos = t.atom_slice(keep).xyz

    centers = physical.avg_pore_loc(4, pos)

    r = physical.limits(pos, centers)[0]

    d = 2*np.mean(r)*math.sin((math.pi / 180) * (args.layer_offset / 2))
    theta = math.atan2(d, args.dbwl) * (180 / math.pi)
    D = math.sqrt(args.dbwl**2 + d**2)

    print 'Angle between stacked rings: %s degrees' % theta
    print 'Distance between stacked benzene centroids: %s nm' % D
    print 'Q space distance: %s A^-1' % (2*math.pi / (10*D))