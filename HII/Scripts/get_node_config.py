#! /usr/bin/env python

import argparse

def initialize():

        parser = argparse.ArgumentParser(description='Determine number of processes for MPI/GPU setup')

        parser.add_argument('-i', '--input', type=str, default='Run.sh')
        args = parser.parse_args()

        return args

if __name__ == "__main__":

	args = initialize()

	with open(args.input, 'r') as f:

		a = []
		for line in f:
			a.append(line)

	node_config = 0

	while a[node_config].count('--tasks-per-node') == 0:
		node_config += 1

	config = a[node_config].split()

	nodes = int(config[2])
	tasks = int(config[3][-1])

	print nodes*tasks
