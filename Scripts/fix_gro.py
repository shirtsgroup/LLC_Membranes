#!/usr/bin/python
# Script to reorder and renumber .gro files 

import argparse

parser = argparse.ArgumentParser(description='Rewrite .gro file so that it matches topology')  # allow input from user

parser.add_argument('-i', '--input', default='HII_packed.gro', help='Name of input file')
parser.add_argument('-s', '--solvent', default='Ar', help='Name of solvent')
parser.add_argument('-o', '--output', default='Ordered.gro', help='Name of reordered output file')

args = parser.parse_args()


def reorder(gro, output):
    a = []
    b = []
    c = []

    count = 1
    hundreds = 0
    for line in gro:
        if line.count('HII') == 1:
            if ((count + 100000*hundreds) % 137) != 0:
                a.append(line[0:5].replace(line[0:5], '{:>5}'.format(str(int((count + 100000*hundreds)/137) + 1))) +
                         line[5:15] + line[15:20].replace(line[15:20], '{:>5}'.format(str(count))) + line[20:len(line)])
                count += 1
            else:
                a.append(line[0:5].replace(line[0:5], '{:>5}'.format(str(int((count + 100000*hundreds)/137))))+line[5:15] +
                         line[15:20].replace(line[15:20], '{:>5}'.format(str(count))) + line[20:len(line)])
                count += 1
        if count >= 100000:
            count = 0
            hundreds += 1

    count1 = (count + 100000*hundreds)/137
    count2 = 0
    for line in gro:
        if line.count('NA') != 0:
            b.append(line[0:5].replace(line[0:5], '{:>5}'.format(count2 + 1 + count1)) + line[5:15] +
                     line[15:20].replace(line[15:20], '{:>5}'.format(str(count + count2))) + line[20:len(line)])
            count2 += 1
        if count >= 100000:
            count = 0
            count2 = 0
            hundreds += 1

    count3 = 0
    for line in gro:
        if line.count('%s' % args.solvent) != 0:
            c.append(line[0:5].replace(line[0:5], '{:>5}'.format(count3 + count2 + count1 + 1)) + line[5:15] +
                     line[15:20].replace(line[15:20], '{:>5}'.format(str(count + count2 + count3))) + line[20:len(line)])
            count3 += 1
        if count >= 100000:
            count = 0
            count2 = 0
            count3 = 0
            hundreds += 1

    f = open(output, 'w')
    f.writelines(['This is a .gro file \n', '%s' % (count + count2 + count3 - 1 + 100000*hundreds) + '\n'])
    for line in a:
        f.writelines([line])
    for line in b:
        f.writelines([line])
    for line in c:
        f.writelines([line])
    f.writelines('%s' % gro[len(gro) - 1])

if __name__ == '__main__':
    f = open(args.input, 'r')

    a = []
    for l in f:
        a.append(l)
    f.close()

    reorder(a, args.output)