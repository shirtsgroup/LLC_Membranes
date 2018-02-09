#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import random
from llclib import transform
from llclib import file_rw
import subprocess


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

    args = parser.parse_args()

    return args


def add_water_placeholder(top):

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


def random_pt_spherical_shell(pt, rmin, rmax):

    # randomly choose a point within the shell created between two spheres of radius rmin and rmax centered at (0,0,0)
    random_pt = np.random.normal(size=3)
    random_pt /= np.linalg.norm(random_pt)
    random_pt *= random.uniform(rmin, rmax)

    # translate to random_pt to be centered at pt
    random_pt += pt

    return random_pt


def random_water_orientation(water_xyz, water_alignment_vector, placement):

    u = np.random.normal(size=3)
    u /= np.linalg.norm(u)

    R = transform.Rvect2vect(water_alignment_vector, u)

    water_xyz -= water_xyz[0, :]  # center at origin

    rotated = np.zeros([water_xyz.shape[1], 3])
    for i in range(water_xyz.shape[1]):
        rotated[i, :] = np.dot(R, water_xyz[i, :])

    rotated += placement

    return rotated


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)
    coordinates = t.xyz[0, :, :]
    natoms = t.n_atoms
    ids = [a.name for a in t.topology.atoms]
    res = [a.residue.name for a in t.topology.atoms]
    full_box = t.unitcell_vectors
    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                    full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]

    water = md.load('/home/bcoscia/PycharmProjects/GitHub/HII/top/solutes/water.gro')
    water_xyz = water.xyz[0, :, :]
    water_centroid = np.mean(water_xyz, axis=0)
    water_alignment_vector = water_xyz[0, :] - water_centroid

    nwater = 1000  # to be determined based on concentration needed

    with open(args.top, 'r') as f:
        top = []
        for line in f:
            top.append(line)

    add_water_placeholder(top)

    ref_atoms = [a.index for a in t.topology.atoms if a.name in args.ref]

    ref_atom_locations = t.xyz[0, ref_atoms, :]

    placement_options = np.linspace(0, len(ref_atoms)-1, len(ref_atoms), dtype=int)

    placement_atom = np.random.choice(placement_options)

    rmin = args.rmin
    rmax = args.rmax

    waters = np.zeros([nwater*3, 3])
    water_ids = ['OW', 'HW1', 'HW2']
    water_res = ['SOL', 'SOL', 'SOL']

    for i in range(nwater):

        placement_atom = np.random.choice(placement_options)  # choose which atom to place water molecule near
        water_placement_point = random_pt_spherical_shell(ref_atom_locations[placement_atom, :], rmin, rmax)  # point near placement atom
        placed_water_coordinates = random_water_orientation(water_xyz, water_alignment_vector, water_placement_point)  # randomly orient water molecule and place at water_placement_point
        coordinates = np.concatenate((coordinates, placed_water_coordinates), axis=0)  # add to full list of coordinates
        ids += water_ids
        res += water_res
        file_rw.write_gro_pos(coordinates, 'water.gro', ids=ids, res=res, box=box_gromacs)
        p1 = subprocess.Popen(["cp", "placeholder_top.top", "top_intermediate.top"])
        p1.wait()
        p2 = subprocess.Popen(["sed", "-i", "-e", "s/PLACEHOLDER/%s/g" % (i + 1), "top_intermediate.top"])
        p2.wait()
        p3 = subprocess.Popen(["gmx", "grompp", "-p", "top_intermediate.top", "-f", "em.mdp", "-o", "em", "-c", "water.gro"])
        p3.wait()
        p4 = subprocess.Popen(["gmx", "mdrun", "-v", "-deffnm", "em"])
        exit()


    ids = ['O', 'H', 'H']*len(ref_atoms)*3
    res = ['HOH', 'HOH', 'HOH']*len(ref_atoms)*3
    file_rw.write_gro_pos(waters, 'water.gro', ids=ids, res=res, box=[8, 8, 8])
    exit()

    with open('test.gro', 'w') as f:

        f.write('HII LLC with water added to tail region\n')
        f.write('%s\n' % natoms)
        gro = header + footer
        for line in gro:
            f.write(line)
