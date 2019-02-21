#!/usr/bin/env python

with open("HX_GMX.top", 'r' ) as f:
	a = []
	for line in f:
		a.append(line)

itp = []

bonds_index = 0
while a[bonds_index].count('[ bonds ]') == 0:
	itp.append(a[bonds_index])
	bonds_index += 1

itp.append(a[bonds_index])
itp.append(a[bonds_index + 1])
bonds_index += 2

while a[bonds_index] != '\n':
	data = a[bonds_index].split()[:3]
	itp.append('{:>6s}{:>6s}{:>6s}\n'.format(data[0], data[1], data[2]))
	bonds_index += 1 

while a[bonds_index].count('[ pairs ]') == 0:
	itp.append(a[bonds_index])
	bonds_index += 1

itp.append(a[bonds_index])
itp.append(a[bonds_index + 1])
bonds_index += 2

while a[bonds_index] != '\n':
	data = a[bonds_index].split()[:3]
	itp.append('{:>6s}{:>6s}{:>6s}\n'.format(data[0], data[1], data[2]))
	bonds_index += 1 

while a[bonds_index].count('[ angles ]') == 0:
	itp.append(a[bonds_index])
	bonds_index += 1

itp.append(a[bonds_index])
itp.append(a[bonds_index + 1])
bonds_index += 2

while a[bonds_index] != '\n':
	data = a[bonds_index].split()[:4]
	itp.append('{:>6s}{:>6s}{:>6s}{:>6s}\n'.format(data[0], data[1], data[2], data[3]))
	bonds_index += 1 

while a[bonds_index].count('[ dihedrals ]') == 0:
	itp.append(a[bonds_index])
	bonds_index += 1

itp.append(a[bonds_index])
itp.append(a[bonds_index + 1])
bonds_index += 2

while a[bonds_index] != '\n':
	data = a[bonds_index].split()[:5]
	itp.append('{:>6s}{:>6s}{:>6s}{:>6s}{:>6s}\n'.format(data[0], data[1], data[2], data[3], data[4]))
	bonds_index += 1 

while a[bonds_index].count('[ dihedrals ]') == 0:
	itp.append(a[bonds_index])
	bonds_index += 1

itp.append(a[bonds_index])
itp.append(a[bonds_index + 1])
bonds_index += 2

while a[bonds_index] != '\n':
	data = a[bonds_index].split()[:5]
	itp.append('{:>6s}{:>6s}{:>6s}{:>6s}{:>6s}\n'.format(data[0], data[1], data[2], data[3], data[4]))
	bonds_index += 1 

with open('test.itp', 'w') as f:
	for line in itp:
		f.write(line)
