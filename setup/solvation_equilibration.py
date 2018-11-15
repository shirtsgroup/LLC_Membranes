#!/usr/bin/env python

import argparse
import subprocess
import sqlite3 as sql
from LLC_Membranes.llclib import topology
import numpy as np
import os

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in


def initialize():

    parser = argparse.ArgumentParser(description='Solvate system and adjust so there is a certain wt % of water')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file to add water to')
    parser.add_argument('-ratio', '--ratio', default=1.5, type=float, help='Ratio of water in pores to water in tails')
    parser.add_argument('-wt', '--weight_percent', default=10, type=float, help='Total weight percent of water')
    parser.add_argument('-tol', '--tolerance', default=1, type=int, help='Number of water molecules')

    # parallelization
    parser.add_argument('-mpi', '--mpi', action="store_true", help='Run MD simulations in parallel')
    parser.add_argument('-np', '--nproc', default=4, help='Number of MPI processes')

    # same flags as to build.py
    parser.add_argument('-b', '--build_monomer', default='NAcarb11V.gro', type=str, help='Name of single monomer'
                        'structure file (.gro format) used to build full system')
    parser.add_argument('-m', '--monomers_per_column', default=20, type=int, help='Number of monomers to stack in each'
                                                                                  'column')
    parser.add_argument('-c', '--ncolumns', default=5, type=int, help='Number of columns used to build each pore')
    parser.add_argument('-r', '--pore_radius', default=.6, type=float, help='Initial guess at pore radius (nm)')
    parser.add_argument('-p', '--p2p', default=4.5, type=float, help='Initial pore-to-pore distance (nm)')
    parser.add_argument('-n', '--nopores', default=4, type=int, help='Number of pores (only works with 4 currently)')
    parser.add_argument('-d', '--dbwl', default=.37, type=float, help='Distance between vertically stacked monomers'
                                                                      '(nm)')
    parser.add_argument('-pd', '--parallel_displaced', default=0, type=float, help='Angle of wedge formed between line'
                        'extending from pore center to monomer and line from pore center to vertically adjacent monomer'
                                                                                   'head group.')
    parser.add_argument('-angles', '--angles', nargs='+', default=[90, 90, 60], type=float, help='Angles between'
                        'box vectors')

    # flags unique to equil.sh
    parser.add_argument('-ring_restraints', '--ring_restraints', default=["C", "C1", "C2", "C3", "C4", "C5"], nargs='+',
                        help='Name of atoms used to restrain head groups during initial equilibration.')

    return parser


if __name__ == "__main__":

    args = initialize().parse_args()

    # get build monomer molecular weight and calculate mw of entire dry system
    build_monomer = topology.Residue(args.build_monomer.split('.')[0])
    if build_monomer.residues:
        build_monomer_mw = sum([x.mw for x in build_monomer.residues])
    else:
        build_monomer_mw = build_monomer.mw

    nmon = args.nopores * args.ncolumns * args.monomers_per_column  # number of build monomers in the system
    dry_mass = nmon * build_monomer_mw  # mass of dry system

    # calculate required water in pores and tails
    water = topology.Residue('HOH')
    args.weight_percent /= 100.0  # convert to fraction
    total_water = int((args.weight_percent * nmon * build_monomer_mw) / (water.mw * (1 - args.weight_percent)))
    tail_water = total_water / (args.ratio + 1)
    pore_water = int(args.ratio * tail_water)
    tail_water = int(tail_water)

    # read database of pore radii and associated water contents
    connection = sql.connect("%s/water_content.db" % location)  # database created in this directory
    crsr = connection.cursor()
    sql_command = "select nwater, radius from radii;"
    sql_output = crsr.execute(sql_command).fetchall()
    nwater = [x[0] for x in sql_output]
    radii = [float(x[1]) / 10 for x in sql_output]
    bin = np.digitize(pore_water, nwater)

    #if bin >= len(nwater):  # will need to implement that

    upper_bound, lower_bound = nwater[bin], nwater[bin - 1]  # upper bound is exclusive. Lower bound is inclusive

    # calculate water content if an entry doesn't already exist
    if abs(pore_water - lower_bound) > args.tolerance:
        # linearly interpolate what the next pore radius should be based on desired amount of water in pore
        interpolation = (pore_water - lower_bound) / (upper_bound - lower_bound)
        r = radii[bin - 1] + (radii[bin] - radii[bin - 1])*interpolation
    else:
        r = radii[bin - 1]

    connection.close()

    subprocess.call(['build.py', '-b', '%s' % args.build_monomer, '-m', '%s' % args.monomers_per_column, '-c', '%s'
                     % args.ncolumns, '-r', '%s' % r, '-p', '%s' % args.p2p, '-n', '%s' % args.nopores,
                     '-d', '%s' % args.dbwl, '-pd', '%s' % args.parallel_displaced])

    subprocess.call(['equil.sh', '-q', '1', '-b', '%s' % args.build_monomer.split('.')[0], '-r',
                     '%s' % ' '.join(args.ring_restraints), '-m', '%s' % args.mpi, '-p', '%s' % args.nproc, ])
    print(nwater)
    exit()