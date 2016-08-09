#!/usr/bin/python
# Script to rewrite and renumber .gro files convert from packmol ouput

import argparse

parser = argparse.ArgumentParser(description = 'Rewrite .gro file so that it matches topology') # allow input from user

parser.add_argument('-i', '--input', default='HII_packed.gro', help = 'Name of input file')

args = parser.parse_args()

no_monomers = 480

f = open(args.input, 'r')

gro = []
for line in f:
	gro.append(line)

a = []
b = []

count = 1
for line in gro:
	if line.count('HII') == 1:
		if (count % 137) != 0:
			a.append(line[0:5].replace(line[0:5], '{:>5}'.format(str(int(count/137) + 1))) + line[5:15] + line[15:20].replace(line[15:20], '{:>5}'.format(str(count))) + line[20:len(line)])
			count += 1
		else:	
                        a.append(line[0:5].replace(line[0:5], '{:>5}'.format(str(int(count/137)))) + line[5:15] + line[15:20].replace(line[15:20], '{:>5}'.format(str(count))) + line[20:len(line)])
                        count += 1

count2 = 1
for line in gro:
	if line.count('NA') != 0:
		b.append(line[0:5].replace(line[0:5], '{:>5}'.format(count2 + count/137)) + line[5:15] + line[15:20].replace(line[15:20], '{:>5}'.format(str(count))) + line[20:len(line)])
		count += 1
		count2 += 1
print len(b)
f = open('packed.gro', 'w')
f.writelines(['This is a .gro file \n', '%s' %(no_monomers*138) + '\n'])
for line in a:
	f.writelines([line])
for line in b:
	f.writelines([line])
f.writelines('%s' %gro[len(gro) - 1])
