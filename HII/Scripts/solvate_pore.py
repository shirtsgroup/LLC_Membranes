#! /usr/bin/env python

import numpy as np
import mdtraj as md
import argparse
from llclib import file_rw
from llclib import physical
import subprocess
import Atom_props
import lc_class


def initialize():

    parser = argparse.ArgumentParser(description='Solvate system and adjust so there is a certain wt % of water')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Initial configuration, preferrably one that has'
                                                                  'been somewhat equilibrate - i.e. the tails wiggled')
    parser.add_argument('-b', '--build_monomer', default='NAcarb11V.gro', help='Name of monomer used to build system')
    parser.add_argument('-w', '--weight_percent', default=6, type=float, help='Desired weight percent of water in pores')
    parser.add_argument('-r', '--radius', default=1, type=float, help='Allowable radius for water to exist from pore center')
    parser.add_argument('-o', '--output', default='solvated.gro', type=str, help='Name of output file')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    subprocess.call(['gmx', 'solvate', '-cs', 'spc216.gro', '-cp', '%s' % args.gro, '-o', 'solv.gro'])

    t = md.load('solv.gro')  # load solvated structure

    mon = lc_class.LC('%s' % args.build_monomer)  # create monomer object
    mw = mon.MW + sum([Atom_props.mass[i] for i in mon.ions])  # molecular weight of monomer with counter ion
    keep = [a.index for a in t.topology.atoms if a.residue.name != 'HOH']

    sol = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']  # only the oxygen atoms
    nsol = len(sol)  # each oxygen atom represents a full SOL residue
    nmon = (t.n_atoms - 3*len(sol)) / (mon.natoms + mon.no_ions)  # number of monomers in the system
    msol = nsol*18.01528  # mass of solutes (molar mass basis)
    mmon = nmon * mw  # mass of monomer (molar mass basis)

    nsol_desired = (mmon * (args.weight_percent/100)) // ((1 - (args.weight_percent/100))*18.01528)

    print('Solvated weight percent: %2.2f %%' % (100*msol/(msol + mmon)))  # fyi

    centers = physical.avg_pore_loc(4, t.xyz[0, keep, :])  # pore centers based on all atoms

    nsol = 0
    pore_solute = []
    for i in sol:
        d = []
        for j in range(centers.shape[1]):
            d.append(np.linalg.norm(centers[:, j] - t.xyz[0, i, :2]))
        if sorted(d)[0] <= args.radius:
            pore_solute.append(i)
            nsol += 1

    msol = nsol*18.01528
    print('Weight percent after restriction to pores: %2.2f %%' % (100*msol/(msol + mmon)))

    remove = int(nsol - nsol_desired)

    if remove < 0:
        print('The pore is not large enough to fit the desired amount of water')
        exit()

    removal_indices = np.zeros([nsol])  # array with number of entries equal to number of remaining solutes
    removal_indices[:remove] = 1  # fill the beginning indices with 1's
    np.random.shuffle(removal_indices)  # shuffle the 1's around to randomize

    for i in range(removal_indices.shape[0]):
        if removal_indices[i] == 0:
            keep.append(pore_solute[i])
            keep.append(pore_solute[i] + 1)
            keep.append(pore_solute[i] + 2)

    final = t.atom_slice(keep)
    nsol = [a.name for a in final.topology.atoms if a.residue.name == 'HOH']
    msol = (len(nsol)/3)*18.01
    print('Final weight percent: %2.2f %%' % (100*msol/(msol + mmon)))  # fyi

    file_rw.write_gro(final, '%s' % args.output)
