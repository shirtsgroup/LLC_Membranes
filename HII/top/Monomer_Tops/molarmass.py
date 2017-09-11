#!/usr/bin/env python

import argparse

def initialize():

    parser = argparse.ArgumentParser(description='Calculate Molar Mass')

    parser.add_argument('-t', '--top', type=str, default='NAcarb11V.itp', help='Topology with masses')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

	args = initialize()

	with open(args.top, 'r') as f:

		a = []
		for line in f:
			a.append(line)

	atoms = 0
	while a[atoms].count('[ atoms ]') == 0:
		atoms += 1

	atoms += 2
	
	mass = []
	while a[atoms] != '\n':
		mass.append(float(a[atoms].split()[-1]))
		atoms += 1

	print('Molar Mass : %s' % sum(mass))
