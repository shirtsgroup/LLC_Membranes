#!/usr/bin/env python

from LLC_Membranes.setup import place_solutes as ps, genmdp
from LLC_Membranes.llclib import topology
import numpy as np
import argparse
import os
import subprocess


def initialize():

    parser = argparse.ArgumentParser(description='Convenience script for adding solutes to pores')

    parser.add_argument('-g', '--gro', default='nosolute.gro', type=str, help='Coordinate file to which solutes will be'
                        'added')
    parser.add_argument('-o', '--out', default='initial.gro', type=str, help='Name of output topology file')
    parser.add_argument('-r', '--solute_resname', default='ETH', type=str, help='Name of solute residue that is being'
                        'added')
    parser.add_argument('-n', '--nsolutes', default=6, type=int, help='Number of solute molecules to add per pore')
    parser.add_argument('-mdps', '--generate_mdps', action="store_true", help='Create input .mdp files')
    parser.add_argument('-noxlink', action="store_false", help='If the system is not cross-linked, add this flag')
    parser.add_argument('-mpi', '--mpi', default=False, help="Specify number of MPI processes")

    return parser


if __name__ == "__main__":

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    args = initialize().parse_args()

    solvent = ps.Solvent('%s' % args.gro, xlink=args.noxlink)
    solute = topology.Solute('%s' % args.solute_resname)

    if args.mpi:
        solvent.mpi = True
        solvent.np = int(args.mpi)

    zbox = solvent.box_vectors[2, 2]
    z = np.linspace(0, zbox, args.nsolutes*2 + 1)[1::2]  # equally space residues
    for i in range(args.nsolutes):
            solvent.place_solute_pores(solute, z=z[i])

    solvent.write_config(name='%s' % args.out)

    if args.generate_mdps:

        mdp = genmdp.SimulationMdp('%s' % args.out, length=5000, barostat='berendsen', xlink=args.noxlink, frames=50)
        mdp.write_em_mdp()
        mdp.write_npt_mdp(out='berendsen')

        mdp = genmdp.SimulationMdp('%s' % args.out, length=1000000, barostat='Parrinello-Rahman', xlink=args.noxlink,
                                   genvel='no', frames=2000)
        mdp.write_npt_mdp(out='PR')

    # put everything in monoclinic cell
    pipe = subprocess.Popen(['echo', '0'], stdout=subprocess.PIPE)
    put_in_box = "gmx trjconv -f %s -o %s -pbc atom -ur tric -s em.tpr" % (args.out, args.out)
    subprocess.Popen(put_in_box.split(), stdin=pipe.stdout)
