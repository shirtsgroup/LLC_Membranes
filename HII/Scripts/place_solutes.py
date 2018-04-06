#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import Atom_props
from llclib import file_rw
import subprocess
import os


def initialize():

    parser = argparse.ArgumentParser(description='Add specified amount of solvent to box')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file to add solutes to')
    parser.add_argument('-c', '--concentration', type=float, help='Concentration of solute (M)')
    parser.add_argument('-n', '--n_solute', type=int, help='Number of solute molecules to add (overrides '
                                                            'concentration')
    parser.add_argument('-cs', '--solute_configuration', help='.gro file for solute molecules')
    parser.add_argument('-ts', '--solute_topology', help='Name of topology file describing solute')

    args = parser.parse_args()

    return args


def random_point_box(box_vectors):
    """

    :param box_vectors: (numpy array, (3, 3)) box vectors. Each row represents a box vector.
    :return: (numpy array, (3)) coordinates of point that lies in box
    """

    A = box_vectors[0, :]  # x box vector
    B = box_vectors[1, :]  # y box vector
    C = box_vectors[2, :]  # z box vector
    u, v, w = np.random.rand(3)  # generate 3 random numbers between 0 and 1
    pt = np.array([0, 0, 0]) + u * A + v * B + w * C  # places point inside 3D box defined by box vector A, B and C

    return pt


def concentration_to_nsolute(conc, box_vectors, solute):
    """
    :param conc: (float) desired solute concentration (M)
    :param box_vectors: (numpy array, (3, 3)) box vectors. Each row represents a box vector.
    :param solute: mdtraj trajectory object generated from solute configuration file (.gro)
    :return: (int) number of solute molecules to add to box to achieve desired concentration
    """

    V = np.dot(box_vectors[2, :], np.cross(box_vectors[0, :], box_vectors[1, :]))  # box volume (nm^3)
    V *= 1 * 10 ** -24  # convert to L
    mols_solute = conc * V  # number of mols of solvent to add

    mw = 0  # molecular weight (grams)
    for a in solute.topology.atoms:
        mw += Atom_props.mass[a.name]

    mass_to_add = mw * mols_solute

    NA = 6.022 * 10 ** 23  # avogadro's number
    mass_solute = mw / NA  # mass of a single solutes (grams)

    nsolute = int(mass_to_add / mass_solute)  # number of solute molecules to add

    actual_concentration = nsolute / (NA*V)  # mol/L

    return nsolute, actual_concentration


def add_water_placeholder(top):
    """
    Open up topology, add water include statement and a PLACEHOLDER so that varying amounts of water can be placed
    :param top: Name of topology to be modified
    :return: placeholder_top.top is written
    """

    solvated = False
    for i in range(len(top)):
        if top[i].count('Water Topology') >= 1:
            solvated = True

    if not solvated:
        # find [ system ] directive
        system_ndx = 0
        while top[system_ndx].count('[ system ]') == 0:
            system_ndx += 1

        top.insert(system_ndx, '; Water Topology\n')
        top.insert(system_ndx + 1, '#include "%s/../top/Forcefields/gaff/tip3p.itp"\n' % location)
        top.insert(system_ndx + 2, '\n')

    top.append('SOL                PLACEHOLDER')

    with open('placeholder_top.top', 'w') as f:

        for line in top:
            f.write(line)


def energy_minimize(steps, nsol, rem, water_placement_point):
    """
    Energy minimize a configuration
    :param steps: number of steepest descent energy minimization steps to take
    :param nsol: number of solute molecules already in the system
    :return: coordinates of energy minimized structure, updated coordinates of reference atoms
    """

    file_rw.write_em_mdp(steps)  # write em.mdp with a given number of steps

    p1 = subprocess.Popen(["cp", "placeholder_top.top", "top_intermediate.top"])  # make a copy of placeholder_top.top
    p1.wait()
    p2 = subprocess.Popen(["sed", "-i", "-e", "s/PLACEHOLDER/%s/g" % nsol, "top_intermediate.top"])  # replace PLACEHOLDER with nwater
    p2.wait()
    p3 = subprocess.Popen(["gmx", "grompp", "-p", "top_intermediate.top", "-f", "em.mdp", "-o", "em", "-c", "water.gro",
                           "-n", "freeze_index.ndx"], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # generate atomic level input file
    p3.wait()
    p4 = subprocess.Popen(["gmx", "mdrun", "-deffnm", "em"], stdout=open(os.devnull, 'w'),
                          stderr=subprocess.STDOUT)  # run energy minimization
    p4.wait()


if __name__ == "__main__":

    args = initialize()

    solvent = md.load(args.gro)
    solvent_box = solvent.unitcell_vectors[0, :, :]

    solute = md.load(args.solute_configuration)

    if args.concentration:
        nsolute, actual_concentration = concentration_to_nsolute(args.concentration, solvent_box, solute)
        print("Actual Concentration : %.2f mol/L" % actual_concentration)
    elif args.n_solute:
        nsolute = args.n_solute
    else:
        print("You must specify a concentration or number of solute molecules")
        exit()

    for n in range(nsolute):
        placement = random_point_box(solvent_box)