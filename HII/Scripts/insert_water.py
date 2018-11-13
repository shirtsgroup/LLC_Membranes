#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from llclib import physical
from llclib import transform
from llclib import file_rw
import subprocess
import lc_class
import top


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Non-solvated initial configuration')
    parser.add_argument('-ref', default='C6', type=str, help='Reference atom for locating pore centers')
    parser.add_argument('-cut', default=0.8, type=float, help='Distance from pore center where water will not be placed')
    parser.add_argument('-pw', '--pore_water', default=4.0, type=float, help='Weight percent of water in pores')
    parser.add_argument('-tw', '--tail_water', default=3.0, type=float, help='Weight percent of water in tail region')
    parser.add_argument('-m', '--monomer', default='NAcarb11V.gro', type=str, help='Monomer system is built with')
    parser.add_argument('-ox', nargs='+', default=['O5', 'O6', 'O7', 'O8', 'O9', 'O10'], help='Oxygen atoms that water'
                                                                                              'will be placed near')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)
    mon = lc_class.LC('NAcarb11V.gro')
    pw = args.pore_water
    tw = args.tail_water

    pore_centers = physical.avg_pore_loc(4, t.xyz)

    # figure out how many waters to each region (pore and tail regions)
    nmon = t.n_atoms / (mon.natoms + mon.no_ions)  # number of monomers in the system
    mon_mass = nmon * mon.MW  # molar mass of unit cell
    mw_water = 18.0154  # molecular weight of water
    sys_mass = mon_mass / (1 - ((pw + tw)/100))  # total molecular weight of system with desired water percentages
    water_mass = sys_mass - mon_mass  # total mass of water in unit cell
    tail_water = int((tw / (pw + tw))*water_mass / mw_water)   # number of waters that belong in the tail region
    pore_water = int((pw / (pw + tw))*water_mass / mw_water)

    pore_no_tails_wt = 100*pore_water*mw_water / (pore_water*mw_water + mon_mass)  # wt % of water in pore region
    # excluding the tail region water. This is necessary to get the proper wt % out of solvate_pore.py

    # subprocess.call(["solvate_pore.py", "-g", "%s" % args.gro, "-b", "%s" % args.monomer, "-w", "%s" % pore_no_tails_wt,
    #                  "-r", "%s" % args.cut, "-o", "solv.gro"])

    t = md.load('solv.gro')

    ids = [a.name for a in t.topology.atoms]
    res_names = [a.residue.name for a in t.topology.atoms]
    # wata = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']
    # pore_water = len(wata)  # The actual number of water molecules placed in the pores after using solvate_pore.py

    placed = t.xyz.shape[1]  # a running total of the number of atoms placed
    positions = np.zeros([placed + 3*tail_water, 3])
    positions[:placed, :] = t.xyz[0, :, :]

    # load up a water molecule
    w = md.load('/home/bcoscia/PycharmProjects/GitHub/HII/top/solutes/water.gro')

    H = w.xyz[0, 1, :]  # location of first hydrogen

    # find places to put water
    ox = [a.index for a in t.topology.atoms if a.name in args.ox]

    # translate a water molecule to a random tail oxygen based on the 1st hydrogen of the water molecule
    positions[placed:(placed + 3), :] = transform.translate(w.xyz[0, :, :], H, t.xyz[0, np.random.choice(ox), :])
    placed += 3
    ids += ["OW", "HW1", "HW2"]
    res_names += 3*["HOH"]

    full_box = t.unitcell_vectors
    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                   full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]

    file_rw.write_gro_pos(positions[:placed, :], 'tail_water.gro', ids=ids, res=res_names, box=box_gromacs)

    # NOTES on what to do next

    # 1) Shift placement of water molecules so that they don't lie directly on top of the randomly chosen oxygen atom
    # 2) Energy minimize a structure with the water placed
    # 3) See how many water molecules can be placed between energy minimizations (all of them?)




