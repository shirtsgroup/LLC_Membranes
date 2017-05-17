#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from llclib import physical
import matplotlib.pyplot as plt
from pymbar import timeseries


def initialize():

    parser = argparse.ArgumentParser(description='Measure the density inside the hydrophobic region of an LLC system')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help = 'Trajectory file (.xtc or .trr). The'
                        'trajectory should be preprocessed beforehand to make molecules whole or else the distance '
                        'calculation get thrown off. Use "gmx trjconv -pbc whole" to do this')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help = 'Name of coordinate file')
    parser.add_argument('-c', '--components', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'],
                        help='component used to track pore positions and define pore radius')
    parser.add_argument('--noshow', action="store_true", help='Do not show the plots')
    parser.add_argument('--save', action="store_true", help='Save plots')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)  # load gromacs trajectory
    keep = [a.index for a in t.topology.atoms if a.name in args.components]   # get the index of all atoms in components
    pore_components = t.atom_slice(keep)  # create a new trajectory object only describing those atoms
    pos = pore_components.xyz

    pcenters = physical.avg_pore_loc(4, pos)  # find the pore centers

    r, r_std = physical.limits(pos, pcenters)

    poresize_equil = timeseries.detectEquilibration(r)[0]
    order_equil = timeseries.detectEquilibration(r/r_std)[0]

    print 'Pore size equilibrated after %d ns' % (t.time[poresize_equil] / 1000)
    print 'Average Pore Size: %.2f +/- %.2f nm' %(np.mean(r[poresize_equil:]), np.mean(r_std[poresize_equil:]))
    print 'Order parameter equilibrated after %d ns' % (t.time[order_equil] / 1000)
    print 'Average Order Parameter: %.2f' % np.mean(r[order_equil:]/r_std[order_equil:])

    plt.figure(1)
    plt.plot(t.time, r)
    plt.title('Pore size vs time')
    plt.ylabel('Average distance of benzene from pore center (nm)')
    plt.xlabel('Time')
    if args.save:
        plt.savefig('poresize.png')

    # plt.figure(2)
    # plt.plot(t.time, r_std)
    # plt.ylabel('Standard deviation of distance of benzene from pore center')
    # plt.xlabel('Time')

    plt.figure(3)
    plt.plot(t.time, r/r_std)
    plt.title('Order parameter vs time')
    plt.xlabel('Time (ps)')
    plt.ylabel('Order parameter')
    if args.save:
        plt.savefig('order.png')

    if not args.noshow:
        plt.show()