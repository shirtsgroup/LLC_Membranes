#usr/bin/env python

import argparse
import numpy as np


def initialize():

    parser = argparse.ArgumentParser(description='Rewrite .gro file so that it matches topology')  # allow input from user

    parser.add_argument('-i', '--input', default='HII_packed.gro', type=str, help='Name of input file')
    parser.add_argument('-s', '--solvent', default='SOL', help='Name of solvent')
    parser.add_argument('-o', '--output', default='Ordered.gro', help='Name of reordered output file')
    parser.add_argument('-n', '--nres', default=480, type=int, help='Number of residues in system')
    parser.add_argumnet('-I', '--ion', default='NA', type=str, help='Name of ion in system')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

	args = intialize()

	with open(args.input, 'r') as f:

		a = []
		for line in f:
			a.append(line)

	atoms_per = (len(a) - 3)/args.nres

	print(atoms_per)
