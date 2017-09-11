#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser(description = 'Run Cylindricity script')

parser.add_argument('-n', '--name', default='Ben')

args = parser.parse_args()

name = args.name

if name[0].isupper():
    if name == 'Ben':
        print 'What a cool guy'
        exit()
    elif name == 'David':
        print 'What a douche'
        exit()
else:
    'Capitalize the first letter dumbass'