#!/usr/bin/env python

with open('NAcarb11V_gromos2016H66.itp', 'r') as f:

	a = []
	for line in f:
		a.append(line)

bonds = 0
while a[bonds].count('[ bonds ]') == 0:
	bonds += 1

bonds += 2

while a[bonds] != '\n':
	
	info = a[bonds].split()
	a[bonds] = '{:>5s}{:>5s}{:>5s}\n'.format(info[0], info[1], info[2])
	bonds += 1

angles = 0
while a[angles].count('[ angles ]') == 0:
        angles += 1

angles += 2

while a[angles] != '\n':

        info = a[angles].split()
        a[angles] = '{:>5s}{:>5s}{:>5s}{:>5s}\n'.format(info[0], info[1], info[2], info[3])
        angles += 1

d = 0
while a[d].count('[ dihedrals ]') == 0:
        d += 1

d += 3

while a[d] != '\n':

        info = a[d].split()
        a[d] = '{:>5s}{:>5s}{:>5s}{:>5s}{:>5s}\n'.format(info[0], info[1], info[2], info[3], info[4])
        d += 1

while a[d].count('[ dihedrals ]') == 0:
        d += 1

d += 2

while a[d] != '\n':

        info = a[d].split()
        a[d] = '{:>5s}{:>5s}{:>5s}{:>5s}{:>5s}\n'.format(info[0], info[1], info[2], info[3], info[4])
        d += 1


with open('new_NAcarb11V_gromos2016H66.itp', 'w') as f:

	for line in a:
		f.write(line)

