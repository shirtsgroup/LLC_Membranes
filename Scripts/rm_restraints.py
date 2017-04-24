#! /usr/bin/env python

import argparse


def initialize():

    parser = argparse.ArgumentParser(description = 'Remove restraints from a topology file')

    parser.add_argument('-p', '--top', default='dipole.itp', type=str, help='Gromacs .itp topology file')
    parser.add_argument('-o', '--out', type=str, help='Output .itp file name. If not specified,'
                                'the output name will be the same as the input name')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    a = []
    with open(args.top) as f:
        for line in f:
            a.append(line)

    count = 0
    while a[count].count('[ position_restraints ]') == 0:
        count += 1

    start = count

    while a[count] != '\n' and count < len(a) - 1:
        count += 1

    end = count + 1

    del a[start:end]

    if args.out:
        with open('%s' % args.out, 'w') as f:
            for line in a:
                f.write(line)
    else:
        with open('%s' % args.top, 'w') as f:
            for line in a:
                f.write(line)

    print 'Restraints removed from topology BIOTCH'