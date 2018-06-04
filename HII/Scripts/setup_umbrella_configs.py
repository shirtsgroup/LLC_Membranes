#!/usr/bin/env python

import numpy as np
import mdtraj as md
import place_solutes
import argparse


def initialize():

    parser = argparse.ArgumentParser(description='Add specified amount of solvent to box')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file to add solutes to')
    parser.add_argument('-n', '--n_configs', default=12, type=int, help='Number of configurations to generate')
    parser.add_argument('-d', '--spacing', type=float, default=0.2, help='Spacing between restraints')
    parser.add_argument('-s', '--solute', default='ETH', help='name of solute residue to add')
    parser.add_argument('-l', '--layers', default=20, type=int, help='number of layers in initial configuration')
    parser.add_argument('-p', '--pores', default=4, type=int, help='number of pores to add solutes to')
    parser.add_argument('-frac', default=0.5, type=float, help='Fraction into membrane to start placing solutes')
    parser.add_argument('-o', '--output', default='solvated', help='Prefix to name of output coordinate files')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    system = place_solutes.Solvent(args.gro)  # system to which to add solute
    zbox = system.box_vectors[2, 2]

    z = np.linspace(zbox * args.frac, zbox*args.frac + args.n_configs*args.spacing - args.spacing, args.n_configs)

    for i, d in enumerate(z):
        if i != 0:  # save almost negligible amount of time..but why waste any?
            system = place_solutes.Solvent(args.gro)  # system to add solute to
        solute = place_solutes.Solute(args.solute)  # object with all relevant solute properties
        system.place_solute_pores(solute, d, layers=args.layers, pores=args.pores)
        system.write_config('%s_%s.gro' % (args.output, i))
