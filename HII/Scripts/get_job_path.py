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

	line = 0
	while a[line].count('StdOut') == 0:
		line += 1

	path = str.strip(a[line].replace('StdOut=', ''))

	print path
