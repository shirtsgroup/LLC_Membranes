#!/usr/bin/env python

# Find the membrane thickness based on a maximum and minimum z position
# Also gives the amount of space needed for water molecules given the amount of space wanted between periodic boundaries
# in the z direction

from __future__ import division
from __future__ import print_function
from past.utils import old_div
import argparse
from llclib import physical
import mdtraj as md
import matplotlib.pyplot as plt
from pymbar import timeseries
import numpy as np


def initialize():

    parser = argparse.ArgumentParser(description = 'Measure Membrane Thickness')

    parser.add_argument('-g', '--gro', type=str, default='wiggle.gro', help = 'Path to input file')
    parser.add_argument('-t', '--traj', type=str, default='traj_whole.xtc', help = 'Trajectory file (.xtc or .trr)')
    parser.add_argument('-w', '--water', default=6, help = 'nm of water between layers')
    parser.add_argument('-m', '--no_monomers', default=6, help= 'number of monomers per layer')
    parser.add_argument('-l', '--layers', default=20, help= 'Number of stacked layers in unit cell')
    parser.add_argument('-p', '--pores', default=4, help= 'Number of pores')
    parser.add_argument('-a', '--atoms', default=137, help='Number of atoms excluding ions')
    parser.add_argument('--trajectory', action="store_true")
    parser.add_argument('--ref_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms which will be used as'
                                                                 'a reference for determining max and min z coordinate')
    parser.add_argument('--save', action="store_true", help='Save plot output')
    parser.add_argument('--noshow', action="store_true", help='Do not display plot')

    args = parser.parse_args()

    return args


def thickness(filename):

    f = open(filename, "r")  # .gro file whose positions of Na ions will be read

    a = []  # list to hold lines of file
    for line in f:
        a.append(line)

    line = 0
    while a[line].count('HII') == 0:
        line += 1

    z = []  # list to hold z positions of all atoms

    benz_carbs = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']

    while a[line].count('HII') != 0:
        if str.strip(a[line][11:15]) in benz_carbs:
            z.append(float(a[line][36:44]))
        line += 1

    z_max = max(z)
    z_min = min(z)
    thick = z_max - z_min

    return thick, z_max, z_min

if __name__ == '__main__':

    args = initialize()
    water_layer = int(args.water)  # nm of water wanted between layers

    if not args.trajectory:
        thick, z_max, z_min = physical.thickness('%s' % args.gro, args.ref_atoms)
        tot_thickness = thick + water_layer
        print('Membrane Thickness: %s nm' % thick)
    else:

        # load trajectory and restrict atoms to reference atoms only
        t = md.load(args.traj, top=args.gro)
        keep = [a.index for a in t.topology.atoms if a.name in args.ref_atoms]
        pos = t.atom_slice(keep).xyz

        # Calculate a trajectory of thicknesses
        thick = physical.thickness('%s' % args.gro, args.ref_atoms, pos)[0]

        # Decide when the system is equilibrated and take the average of thicknesses after equilibration
        equil_frame = timeseries.detectEquilibration(thick)[0]
        print('Equilibration detected after %d ns' % (old_div(t.time[equil_frame], 1000.0)))
        print('Average membrane thickness: %.2f +/- %.2f nm' % (np.mean(thick[equil_frame:]),
                                                                  np.std(thick[equil_frame:])))

        # Plot thickness vs time
        plt.plot(t.time, thick)
        plt.title('Membrane thickness vs. time')
        plt.xlabel('Time (ps)')
        plt.ylabel('Membrane thickness (nm)')

        # Save plot if you want
        if args.save:
            plt.savefig('thickness.png')

        # Don't show plot if you don't want to
        if not args.noshow:
            plt.show()
