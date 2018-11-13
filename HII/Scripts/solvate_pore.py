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
    parser.add_argument('-r', '--radius', default=10, type=float, help='Allowable radius for water to exist from pore center (A)')
    parser.add_argument('-o', '--output', default='solvated.gro', type=str, help='Name of output file')
    parser.add_argument('-l', '--gap', type=float, help='Width of gap (if one exists). If there is a gap, ignore those '
                                                        'waters for wt % calculations')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    subprocess.call(['gmx', 'solvate', '-cs', 'spc216.gro', '-cp', '%s' % args.gro, '-o', 'solv.gro'])
    echo = subprocess.Popen(['echo', '0'], stdout=subprocess.PIPE)
    subprocess.Popen(['gmx', 'trjconv', '-f', 'solv.gro', '-o', 'solv.gro', '-s', 'solv.gro', '-ur', 'tric', '-pbc',
                     'atom'], stdin=echo.stdout)
    r = args.radius / 10

    t = md.load('solv.gro')  # load solvated structure
    z = t.unitcell_vectors[0, 2, 2]

    mon = lc_class.LC('%s' % args.build_monomer)  # create monomer object
    mw = mon.MW + sum([Atom_props.mass[i] for i in mon.ions])  # molecular weight of monomer with counter ion
    keep = [a.index for a in t.topology.atoms if a.residue.name != 'HOH']

    if args.gap:
        upper_limit = z - (args.gap / 2)
        lower_limit = args.gap / 2
        sol = []
        sol_ignore = []
        for a in t.topology.atoms:
            if a.residue.name == 'HOH' and a.name == 'O':
                if upper_limit > t.xyz[0, a.index, 2] > lower_limit:
                    sol.append(a.index)
                else:
                    sol_ignore.append(a.index)
                    sol_ignore.append(a.index + 1)
                    sol_ignore.append(a.index + 2)
        # sol = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O' \
        #        and upper_limit > t.xyz[0, a.index, 2] > lower_limit]  # only the oxygen atoms within membrane
    else:
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
        if sorted(d)[0] <= r:
            pore_solute.append(i)
            nsol += 1

    msol = nsol*18.01528
    print('Weight percent after restriction to pores: %2.2f %%' % (100*msol/(msol + mmon)))

    remove = int(nsol - nsol_desired)
    removal_indices = np.zeros([nsol])  # array with number of entries equal to number of remaining solutes

    if remove < 0:
        print('The pore is not large enough to fit the desired amount of water..proceeding anyways')
    else:
        removal_indices[:remove] = 1  # fill the beginning indices with 1's
        np.random.shuffle(removal_indices)  # shuffle the 1's around to randomize

    if args.gap:
        keep += sol_ignore

    nsol = 0
    for i in range(removal_indices.shape[0]):
        if removal_indices[i] == 0:
            nsol += 1
            keep.append(pore_solute[i])
            keep.append(pore_solute[i] + 1)
            keep.append(pore_solute[i] + 2)

    res = []
    ids = []
    for a in t.topology.atoms:
        if a.residue.name == 'HOH':
            res.append('SOL')
            if a.name == 'H1':
                ids.append('HW1')
            elif a.name == 'H2':
                ids.append('HW2')
            elif a.name == 'O':
                ids.append('OW')
        else:
            res.append(a.residue.name)
            ids.append(a.name)

    res = np.array(res)
    ids = np.array(ids)
    full_box = t.unitcell_vectors
    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]  # convert to .gro format

    final_coords = t.xyz[0, keep, :]
    # final = t.atom_slice(keep)
    # nsol = [a.name for a in final.topology.atoms if a.residue.name == 'HOH']
    # msol = (len(nsol)/3)*18.01
    msol = nsol*18.01
    print('Final weight percent: %2.2f %%' % (100*msol/(msol + mmon)))  # fyi
    # print('^^^^ Without water in the tails ^^^^')
    # print('n_water_tails | wt % water pores | wt % water tails | total wt % water')
    # tail_water = [i * 200 for i in range(13)]
    # for w in tail_water:
    #     mtails = w * 18.01
    #     mtot = mtails + msol + mmon
    #     percent_pores = 100 * msol / mtot
    #     percent_tails = 100 * mtails / mtot
    #     percent_tot = 100*(msol + mtails) / mtot
    #     print('{:^14d}|{:^18.2f}|{:^18.2f}|{:^15.2f}'.format(w, percent_pores, percent_tails, percent_tot))

    file_rw.write_gro_pos(final_coords, args.output, box=box_gromacs, res=list(res[keep]), ids=list(ids[keep]))
