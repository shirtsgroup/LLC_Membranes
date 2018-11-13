#!/usr/bin/env python

import argparse
import os
from LLC_Membranes.llclib import topology
import numpy as np
import subprocess
import os

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def initialize():

    parser = argparse.ArgumentParser(description='Balance charge for a given residue by rescaling')

    parser.add_argument('-i', '--itp', help='Name of .itp to be modified (without extension)')
    parser.add_argument('-nc', '--net_charge', default=0, type=int, help='Desired net charge on molecule')
    parser.add_argument('-p', '--precision', default=6, type=int, help='Number of decimal points')

    return parser


if __name__ == "__main__":

    args = initialize().parse_args()

    res = topology.Residue(args.itp)

    charges = res.charges

    for key, value in charges.items():
        charges[key] = int(value * 10**args.precision)

    net_charge = sum(res.charges.values())

    if net_charge == args.net_charge:
        print('Charge is already balanced')
        exit()

    # increment all charge values if net_charge is greater than res.natoms
    n_changes = abs(int(net_charge / res.natoms))

    if net_charge < args.net_charge:
        increment = 1
    else:
        increment = -1

    for i in range(n_changes):
        for key in charges.keys():
            charges[key] += increment
        net_charge += increment*res.natoms

    if net_charge < args.net_charge:
        changes = np.random.choice(list(charges.keys()), size=abs(net_charge - args.net_charge), replace=False)
        for i in changes:
            charges[i] += increment

    subprocess.Popen(['cp', '%s/../top/topologies/%s.itp' % (script_location, args.itp),
                      '%s/../top/topologies/%s.itp.bak' % (script_location, args.itp)])

    with open('%s/../top/topologies/%s.itp' % (script_location, args.itp), 'w') as f:
        line = 0
        while res.itp[line].count('[ atoms ]') == 0:
            f.write(res.itp[line])
            line += 1
        f.write(res.itp[line])
        f.write(res.itp[line + 1])
        line += 2
        while res.itp[line] != '\n':
            data = res.itp[line].split()
            f.write('{:>6d}{:>5s}{:>6d}{:>6s}{:>6s}{:>5d}{:>13.6f}{:>13.6f}\n'.format(int(data[0]), data[1], int(data[2]),
                    data[3], data[4], int(data[5]), charges[data[4]] / 10**args.precision, float(data[7])))
            line += 1

        while line < len(res.itp):
            f.write(res.itp[line])
            line += 1
