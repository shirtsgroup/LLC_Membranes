#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='wiggle.trr', help = 'Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help = 'Some kind of configuration file that mdtraj needs'
                                                                    'Can be a .pdb as well')
    parser.add_argument('--single_frame', help='Specify this flag if you want to analyze one frame only',
                        action="store_true")

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    if args.single_frame:
        t = md.load('%s' % args.gro)
    else:
        t = md.load('%s' % args.traj, top='%s' % args.gro)
    print t

    hbonds = md.wernet_nilsson(t, exclude_water=False)
    print len(hbonds)
    print hbonds[0]
    exit()

    label = lambda hbond : '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))

    for hbond in hbonds:
        print label(hbond)