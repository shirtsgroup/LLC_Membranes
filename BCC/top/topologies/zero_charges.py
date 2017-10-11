#!/usr/bin/env python

with open('BCC_cleaned.itp', 'r') as f:

	a = []
	for line in f:
		a.append(line)

atoms_section = 0
while a[atoms_section].count('[ atoms ]') == 0:
	atoms_section += 1

atoms_section += 2
while a[atoms_section] != '\n':
	a[atoms_section] = a[atoms_section].replace(str.strip(a[atoms_section][38:47]), '0.00000')
	atoms_section += 1

with open('BCC_nocharge.itp', 'w') as f:
	for line in a:
		f.write(line)
