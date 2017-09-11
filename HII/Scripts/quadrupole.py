#! /usr/bin/env python

import numpy as np
import mdtraj as md
import argparse
import lc_class


def initialize():

    parser = argparse.ArgumentParser(description='Write .mdp files and topology files for various types of simulation')

    parser.add_argument('-g', '--gro', default='initial.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-b', '--build_mon', default='NAcarb11V.gro', type=str, help='Name of coordinate file used to'
                                                                                     'build assembly')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    properties = lc_class.LC('%s' % args.build_mon)

    t = md.load('%s' % args.gro)