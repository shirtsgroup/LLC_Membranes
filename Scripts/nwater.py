#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help = 'Configuration file, .pdb or .gro')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)

    water = [a.index for a in t.topology.atoms if a.name == 'O' and 'HOH' in str(a.residue)]

    print len(water)