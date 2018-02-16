#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import random
from llclib import transform
from llclib import file_rw
import subprocess
import copy
import os, glob
import sys
import time


def initialize():

    parser = argparse.ArgumentParser(description='Add water to the tail region of a HII LLC structure')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of structure file where water will be added')
    parser.add_argument('-p', '--top', default='topol.top', help='Name of topology that needs to be updated')
    parser.add_argument('-o', '--out', default='solv_tails.gro', help='Name of output file')
    parser.add_argument('-r', '--ref', nargs='+', default=['O5', 'O6', 'O7', 'O8', 'O9', 'O10'], help='Names of atoms'
                        'where water molecules will be placed in proximity to')
    parser.add_argument('-rmin', default=0.2, type=float, help='Minimum distance away from reference atom to place'
                                                               'water molecules')
    parser.add_argument('-rmax', default=0.4, type=float, help='Maximum distance away from reference atom to place'
                                                               'water molecules')
    parser.add_argument('-n', '--nwater', default=600, type=int, help='Number of waters to add to the system')

    args = parser.parse_args()

    return args


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
        top.insert(system_ndx + 1, '#include "/home/bcoscia/PycharmProjects/GitHub/HII/top/Forcefields/gaff/tip3p.itp"\n')
        top.insert(system_ndx + 2, '\n')

    top.append('SOL                PLACEHOLDER')

    with open('placeholder_top.top', 'w') as f:

        for line in top:
            f.write(line)


def write_em_mdp(steps):
    """
    Write energy minimization .mdp file
    :param steps: number of steps to take using steepest descents algorithm
    :return: Directly writes an energy minimization .mdp file
    """

    with open('em.mdp', 'w') as f:

        f.write("title = Energy Minimization\n")
        f.write("integrator = steep\n")
        f.write("nsteps = %s\n" % steps)
        f.write("cutoff-scheme = verlet\n")
        f.write("nstlist = 40\n")


def random_pt_spherical_shell(pt, rmin, rmax):
    """
    Pick a random point in the region bounded by two spheres
    :param pt: center of concentric spheres
    :param rmin: inside sphere radius
    :param rmax: outside sphere radius
    :return: randomly chosen point between spheres with radii rmin and rmax
    """

    # randomly choose a point within the shell created between two spheres of radius rmin and rmax centered at (0,0,0)
    random_pt = np.random.normal(size=3)  # random vector chosen from gaussian since gaussian is spherically symmetric
    random_pt /= np.linalg.norm(random_pt)  # normalize
    random_pt *= random.uniform(rmin, rmax)  # randomly choose r distance between rmin and rmax

    # translate to random_pt to be centered at pt
    random_pt += pt

    return random_pt


def random_water_orientation(water_xyz, water_alignment_vector, placement):
    """
    Randomly orient a water molecule and then place it a desired location
    :param water_xyz: 3d coordinates of a water molecule
    :param water_alignment_vector: A reference vector to rotate the water molecule about
    :param placement: where to place final water configuration in space
    :return: coordinates of oriented and translated water molecule
    """

    u = np.random.normal(size=3)  # random vector. From normal distribution since sphere
    u /= np.linalg.norm(u)  # normalize

    R = transform.Rvect2vect(water_alignment_vector, u)  # rotation matrix to align water_alignment_vector with u

    water_xyz -= water_xyz[0, :]  # center at origin

    rotated = np.zeros([water_xyz.shape[1], 3])
    for i in range(water_xyz.shape[1]):
        rotated[i, :] = np.dot(R, water_xyz[i, :])

    rotated += placement  # translate to deisred location

    return rotated


def energy_minimize(steps, nwater):
    """
    Energy minimize a configuration
    :param steps: number of steepest descent energy minimization steps to take
    :param nwater: number of water molecules in the system
    :return: coordinates of energy minimized structure, updated coordinates of reference atoms
    """

    write_em_mdp(steps)  # write em.mdp with a given number of steps

    p1 = subprocess.Popen(["cp", "placeholder_top.top", "top_intermediate.top"])  # make a copy of placeholder_top.top
    p1.wait()
    p2 = subprocess.Popen(["sed", "-i", "-e", "s/PLACEHOLDER/%s/g" % (nwater), "top_intermediate.top"])  # replace PLACEHOLDER with nwater
    p2.wait()
    p3 = subprocess.Popen(["gmx", "grompp", "-p", "top_intermediate.top", "-f", "em.mdp", "-o", "em", "-c",
                           "water.gro"], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # generate atomic level input file
    p3.wait()
    p4 = subprocess.Popen(["gmx", "mdrun", "-deffnm", "em"], stdout=open(os.devnull, 'w'),
                          stderr=subprocess.STDOUT)  # run energy minimization
    p4.wait()

    t = md.load('em.gro')
    minimized_coordinates = t.xyz[0, :, :]  # coordinates of energy minimized system
    new_ref_atom_locations = t.xyz[0, ref_atoms, :]  # new coordinates of reference atoms

    return minimized_coordinates, new_ref_atom_locations


def place_water(placement_options, ref_atom_locations, coordinates, rmin, rmax, ids, res, nwater, steps):
    """
    Place water molecules near a random reference atom and perform a short energy minimization
    :param placement_options: All of the possible indexes of reference atoms near which a water molecule could be placed (list)
    :param ref_atom_locations: xyz coordinates of reference atoms (numpy array [number of ref atoms, 3]
    :param coordinates: Coordinates of full system (numpy array [natoms, 3])
    :param rmin: minimimum distance to place water molecules from reference atom (float)
    :param rmax: maximum distance to place water molecules from reference atom (float)
    :param ids: all atom names in order they appear in coordinates (list)
    :param res: all residue atoms in order they appear in coordinates (list)
    :param nwater: number of water molecules in the system (int)
    :param steps: number of energy minimization steps to take (int)
    :return: Potential energy of slightly minimized system, coordinates of minimized system, atom used for placement,
    new locations of reference atoms
    """

    placement_atom = np.random.choice(placement_options)  # choose which atom to place water molecule near
    water_placement_point = random_pt_spherical_shell(ref_atom_locations[placement_atom, :], rmin, rmax)  # point near placement atom
    placed_water_coordinates = random_water_orientation(water_xyz, water_alignment_vector, water_placement_point)  # randomly orient water molecule and place at water_placement_point
    new_coordinates = np.concatenate((coordinates, placed_water_coordinates), axis=0)  # add to full list of coordinates
    names = ids + water_ids  # add water to ids
    residues = res + water_res  # add water residue to res
    file_rw.write_gro_pos(new_coordinates, 'water.gro', ids=names, res=residues, box=box_gromacs)  # write out config with new water
    minimized_coordinates, new_ref_atom_locations = energy_minimize(steps, nwater + 1)  # energy minimzed system
    nrg = subprocess.check_output(["awk", "/Potential Energy/ {print $4}", "em.log"])  # get Potential energy from em.log

    return float(nrg.decode("utf-8")), minimized_coordinates, placement_atom, new_ref_atom_locations


if __name__ == "__main__":

    args = initialize()

    if os.path.exists(args.out):  # Get rid of output file if it already exists since it can mess things up
        os.remove(args.out)

    t = md.load(args.gro)
    coordinates = t.xyz[0, :, :]  # initial coordinates
    natoms = t.n_atoms  # number of atoms in the system

    # handle mdtraj renaming SOL to HOH
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

    full_box = t.unitcell_vectors   # unitcell vectors
    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                    full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]  # convert to .gro format

    water = md.load('/home/bcoscia/PycharmProjects/GitHub/HII/top/solutes/water.gro')  # load water structure
    water_xyz = water.xyz[0, :, :]  # get water coordinates
    water_centroid = np.mean(water_xyz, axis=0)  # get centroid of water
    water_alignment_vector = water_xyz[0, :] - water_centroid  # vector around which water molecule can be rotated

    nwater = args.nwater  # to be determined based on concentration needed

    with open(args.top, 'r') as f:
        top = []
        for line in f:
            top.append(line)

    add_water_placeholder(top)  # add a placeholder field to topol.top

    ref_atoms = [a.index for a in t.topology.atoms if a.name in args.ref]  # indices of reference atoms

    ref_atom_locations = t.xyz[0, ref_atoms, :]  # coordinates of reference atoms

    placement_options = [i for i in range(len(ref_atoms))]  # list of indices corresponding to atoms in ref_atom_locations

    rmin = args.rmin  # min distance from ref atom to place water molecule
    rmax = args.rmax  # max distance from ref atom to place water molecule

    water_ids = ['OW', 'HW1', 'HW2']  # water atom names
    water_res = ['SOL', 'SOL', 'SOL']  # water residue name

    write_em_mdp(5)
    trials = 0  # number of water placement tries

    start = time.time()  # when water placement starts

    for i in range(nwater):
        # randomly place water molecule and do short energy minimization
        energy, new_coordinates, placement_atom, new_ref_atom_locations = place_water(placement_options,
                                                        ref_atom_locations, coordinates, rmin,  rmax, ids, res, i, 5)
        trials += 1
        count = 0
        while energy > -50000:  # make sure the potential energy doesn't get too close to exploding

            energy, new_coordinates, placement_atom, new_ref_atom_locations = place_water(placement_options,
                                                            ref_atom_locations, coordinates, rmin, rmax, ids, res, i, 5)
            trials += 1
            count += 1
            if count > 10:
                # if the energy is too high for water placement, run a longer energy minimization
                sys.stdout.write("\r Running longer energy minimization...                                             "
                                 "             \r")
                sys.stdout.flush()
                file_rw.write_gro_pos(coordinates, 'water.gro', box=box_gromacs, ids=ids, res=res)
                coordinates, ref_atom_locations = energy_minimize(500, i)  # energy minimize and update coordinates
                count = 0
                write_em_mdp(1)

        success_rate = 100*((i + 1) / trials)
        s = "Waters placed: %s/%s, Placement Success Rate : %3.2f, Potential Energy = %s\r" % (i, nwater, success_rate, energy)
        sys.stdout.write("\r"+s)
        sys.stdout.flush()
        ids += water_ids
        res += water_res
        coordinates = copy.deepcopy(new_coordinates)
        placement_options.remove(placement_atom)
        for filename in glob.glob("./#*"):
            os.remove(filename)

    print("Waters placed in %4.2f seconds" % (time.time() - start))
    file_rw.write_gro_pos(coordinates, args.out, box=box_gromacs, res=res, ids=ids)