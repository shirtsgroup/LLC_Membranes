#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import top
import lc_class
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Remove bonds/angles/torsions from crosslinks across z PBCs')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of gro whose coordinates we are interested in.'
                        'Leave dummy atoms in the .gro. Preprocess it with gmx trjconv -pbc nojump using a reference'
                        'structure. A good reference structure is the initial structure before crosslinking')
    parser.add_argument('-t', '--top', default='topol.top', help='name of topology')
    parser.add_argument('-b', '--build_mon', default='NAcarb11V.gro', help='Name of monomer system was built with')
    parser.add_argument('-o', '--output', default='capped.itp', help='Name of output topology')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    traj = md.load(args.gro)  # load positions
    pos = traj.xyz
    t = top.Top(args.top, vsites=False)  # load topology
    mon = lc_class.LC(args.build_mon)

    c1 = [a.index + 1 for a in traj.topology.atoms if a.name in mon.c1_atoms]
    c2 = [a.index + 1 for a in traj.topology.atoms if a.name in mon.c2_atoms]

    xlinks = t.get_bonds_between(c1, c2)

    cross_z = []
    for i in range(len(xlinks)):
        d = abs(pos[0, xlinks[i][0], 2] - pos[0, xlinks[i][1], 2])
        if d > 2:
            cross_z.append(xlinks[i][0])
            cross_z.append(xlinks[i][1])

    bonds = t.get_bonds_to([30772, 28769])

    for i in bonds.keys():
        for j in bonds[i]:
            print(j, t.get_type(j))

    exit()
    t.change_atom_type(cross_z, 'cg')

    t.write_top('%s' % args.output)

    #carbons = mon.c1_atoms + mon.c2_atoms  # carbons involved in crosslinking for this system

    # find indices of c1_atoms and c2_atoms separately. Make a function in top.py to do it. Find out bonds

    # topology = []
    # with open(args.top, 'r') as f:
    #     for line in f:
    #         topology.append(line)


