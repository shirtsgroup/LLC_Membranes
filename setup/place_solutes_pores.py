#!/usr/bin/env python

from LLC_Membranes.setup import place_solutes as ps, genmdp
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
    parser.add_argument('-mdps', '--generate_mdps', action="store_true", help='Create input .mdp files')
    parser.add_argument('-noxlink', action="store_false", help='If the system is not cross-linked, add this flag')
    parser.add_argument('-mpi', '--mpi', default=False, help="Specify number of MPI processes")

    return parser


if __name__ == "__main__":

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    args = initialize().parse_args()

    solvent = ps.Solvent('%s' % args.gro, xlink=args.noxlink)
    solute = ps.Solute('%s' % args.solute_resname)

    if args.mpi:
        solvent.mpi = True
        solvent.np = int(args.mpi)

    zbox = solvent.box_vectors[2, 2]
    z = np.linspace(0, zbox, args.nsolutes*2 + 1)[1::2]  # equally space residues
    for i in range(args.nsolutes):
            solvent.place_solute_pores(solute, z=z[i])

    solvent.write_config(name='%s' % args.out)

    if args.generate_mdps:

        mdp = genmdp.SimulationMdp('%s' % args.out, length=5000, barostat='berendsen', xlink=args.noxlink)
        mdp.write_em_mdp()
        mdp.write_npt_mdp(out='berendsen')

        mdp = genmdp.SimulationMdp('%s' % args.out, length=400000, barostat='Parrinello-Rahman', xlink=args.noxlink,
                                   genvel='no')
        mdp.write_npt_mdp(out='PR')

