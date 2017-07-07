#! /usr/bin/env python

import argparse

def initialize():

	parser = argparse.ArgumentParser(description='Get job ids from the queue')

	parser.add_argument('-i', '--input', type=str, default='jobs.py')
	args = parser.parse_args()

	return args

if __name__ == "__main__":

	args = initialize()

	with open(args.input, 'r') as f:
		a = []
		for line in f:
			a.append(line)

	ids = []
	for i in range(1, len(a)):

		ids.append(a[i].split()[0])

	for i in ids:
		print int(i)
