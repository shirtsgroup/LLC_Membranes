#!/usr/bin/env python

import argparse

def initialize():
	
	parser = argparse.ArgumentParser(description='format v3d for use with XRDplot')

	parser.add_argument('-i', '--input', default='output.v3d', type=str, help='Name of input .v3d file')
	parser.add_argument('-o', '--output', default='output.v3d', type=str, help='Name of newly formatted .v3d file')
	parser.add_argument('-x', '--lx', type=float, help='x dimension')
	parser.add_argument('-y', '--ly', type=float, help='y dimension')
	parser.add_argument('-z', '--lz', type=float, help='z dimension')

	args = parser.parse_args()

	return args

if __name__ == "__main__":

	args = initialize()

	with open(args.input, 'r') as f:
		Nx, Ny, Nz = [int(x) for x in next(f).split()]
		I = [float(x) for x in next(f).split()]

	with open(args.output, 'w') as f:
		f.write('{:<4d}{:<4d}{:<4d}\n'.format(Nx, Ny, Nz))
		f.write('{:<3.5f}{:1}{:<3.5f}{:1}{:<3.5f}\n'.format(args.lx * 10,'', args.ly * 10,'', args.lz * 10))
		for intensity in I:
			f.write(str(intensity) + '\n')
	



