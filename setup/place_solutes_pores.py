#!/usr/bin/env python

from LLC_Membranes.setup import place_solutes as ps
import numpy as np
import argparse
import os


def initialize():

    parser = argparse.ArgumentParser(description='Convenience script for adding solutes to pores')

    parser.add_argument('-g', '--gro', default='nosolute.gro', type=str, help='Coordinate file to which solutes will be'
                        'added')
    parser.add_argument('-o', '--out', default='initial.gro', type=str, help='Name of output topology file')
    parser.add_argument('-r', '--solute_resname', default='ETH', type=str, help='Name of solute residue that is being'
                        'added')
    parser.add_argument('-n', '--nsolutes', default=6, type=int, help='Number of solute molecules to add')

    return parser


if __name__ == "__main__":

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    args = initialize().parse_args()

    solvent = ps.Solvent('%s' % args.gro)
    solute = ps.Solute('%s' % args.solute_resname)

    zbox = solvent.box_vectors[2, 2]
    z = np.linspace(0, zbox, args.nsolutes*2 + 1)[1::2]  # equally space residues

    for i in range(args.nsolutes):
            solvent.place_solute_pores(solute, z=z[i])

    solvent.write_config(name='%s' % args.out)
