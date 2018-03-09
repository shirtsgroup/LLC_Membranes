#! /usr/bin/env python

import argparse
import mdtraj as md
import os
import numpy as np
from llclib import physical
from llclib import file_rw
import solvate_tails
import subprocess
from scipy import spatial

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in


def initialize():

    parser = argparse.ArgumentParser(description='Solvate system and adjust so there is a certain wt % of water')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file to add water to')
    parser.add_argument('-t', '--top', default='topol.top', help='Name of topology file to be modified')
    parser.add_argument('-pa', '--pore_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='atoms that'
                        'will be used to define the pore center')
    parser.add_argument('-ta', '--tail_atoms', nargs='+', default=['O5', 'O6', 'O7', 'O8', 'O9', 'O10'], help='Names of'
                        'atoms where water molecules will be placed in proximity to')
    parser.add_argument('-p', '--probability', type=float, default=0.6, help='Probability of water placement in pores. '
                        'The probability to place a water in the tails is therefore 1 - args.probability')
    parser.add_argument('-rp', '--radius_pore', default=0.8, type=float, help='Allowable radius for water to exist from'
                        ' pore center (nm)')
    parser.add_argument('--emsteps', type=int, default=5, help='Number of energy minimization steps to take after'
                        'water insertion')
    parser.add_argument('-rmin', default=0.2, type=float, help='Minimum distance away from reference atom to place'
                        'water molecules')
    parser.add_argument('-rmax', default=0.4, type=float, help='Maximum distance away from reference atom to place'
                        'water molecules')
    parser.add_argument('-rem', default=1, type=float, help='Atoms within this radius of placed water molecules will'
                        'be energy minimized. All others will be frozen (nm)')

    args = parser.parse_args()

    return args


def random_pt_cylinder(R, H):
    """
    Generate a random [x, y, z] point within a cylinder aligned with the z-axis
    :param r: cylinder radius (nm)
    :param h: cylinder height (nm)
    :return: random point in cylinder
    """

    s = np.random.rand()  # random value from uniform distribution over [0, 1)
    theta = 2*np.pi*np.random.rand()  # random value from uniform distribution over [0, 2pi)
    r = np.sqrt(s)*R
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = H*np.random.rand()  # random value from uniform distribution over [0, h)

    return np.array([x, y, z])


def get_ids_res(t):

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

    return ids, res


def water_parameters(water):

    water_xyz = water.xyz[0, :, :]  # get water coordinates
    water_centroid = np.mean(water_xyz, axis=0)  # get centroid of water
    water_alignment_vector = water_xyz[0, :] - water_centroid  # vector around which water molecule can be rotated

    water_ids, water_res = get_ids_res(water)

    return water_xyz, water_centroid, water_alignment_vector, water_ids, water_res


def freeze(water_placement_point, rem):
    """
    Write an index file for atoms to be frozen
    :param water_placement_point: xyz position of where water molecule was placed
    :param rem: spherical radius measured from water molecule placement point outside which all atoms will be frozen
    :return: index file with indices of atoms to be frozen
    """

    t = md.load('water.gro')
    pos = t.xyz

    pts = spatial.cKDTree(pos[0, :, :]).query_ball_point(water_placement_point, rem)

    freeze_indices = [a.index for a in t.topology.atoms if a.index not in pts]

    with open('freeze_index.ndx', 'w') as f:

        f.write('[ Freeze ]\n')
        for i, entry in enumerate(freeze_indices):
            if (i + 1) % 15 == 0:
                f.write('{:5d}\n'.format(entry + 1))
            else:
                f.write('{:5d} '.format(entry + 1))


def energy_minimize(steps, nwater, rem, water_placement_point):
    """
    Energy minimize a configuration
    :param steps: number of steepest descent energy minimization steps to take
    :param nwater: number of water molecules in the system
    :param rem: spherical radius measured from water molecule place point inside which all atoms will be energy minimized
    :param water_placement_point: xyz coordinates where water molecule centroid was placed
    :return: coordinates of energy minimized structure, updated coordinates of reference atoms
    """

    freeze(water_placement_point, rem)  # write index file with group specifiying atoms to be frozen

    file_rw.write_em_mdp(steps, freeze=True, freeze_group='Freeze', freeze_dim='xyz')  # write em.mdp with a given number of steps. Should move this function to file_rw

    p1 = subprocess.Popen(["cp", "placeholder_top.top", "top_intermediate.top"])  # make a copy of placeholder_top.top
    p1.wait()
    p2 = subprocess.Popen(["sed", "-i", "-e", "s/PLACEHOLDER/%s/g" % (nwater), "top_intermediate.top"])  # replace PLACEHOLDER with nwater
    p2.wait()
    p3 = subprocess.Popen(["gmx", "grompp", "-p", "top_intermediate.top", "-f", "em.mdp", "-o", "em", "-c", "water.gro",
                           "-n", "freeze_index.ndx"], stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # generate atomic level input file
    p3.wait()
    p4 = subprocess.Popen(["gmx", "mdrun", "-deffnm", "em"], stdout=open(os.devnull, 'w'),
                          stderr=subprocess.STDOUT)  # run energy minimization
    p4.wait()


def place_water(xyz, water_placement_point, ids, res, box, emsteps, rem):

    water = md.load('%s/../top/solutes/water.gro' % location)  # load water structure
    water_xyz, water_centroid, water_alignment_vector, water_ids, water_res = water_parameters(water)

    placed_water_coordinates = solvate_tails.random_water_orientation(water_xyz, water_alignment_vector,
                                                                      water_placement_point)
    new_coordinates = np.concatenate((xyz, placed_water_coordinates), axis=0)  # add to full list of coordinates
    names = ids + water_ids  # add water to ids
    residues = res + water_res  # add water residue to res
    file_rw.write_gro_pos(new_coordinates, 'water.gro', ids=names, res=residues, box=box)  # write out config with new water
    energy_minimize(emsteps, nwater + 1, rem, water_placement_point)  # energy minimzed system
    nrg = subprocess.check_output(["awk", "/Potential Energy/ {print $4}", "em.log"])  # get Potential energy from em.log

    return float(nrg.decode("utf-8"))

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)
    ids, res = get_ids_res(t)  # get ids and residues
    nwater = res.count('SOL') // 3

    box = t.unitcell_vectors   # unitcell vectors
    box_gromacs = [box[0, 0, 0], box[0, 1, 1], box[0, 2, 2], box[0, 0, 1], box[0, 2, 0], box[0, 1, 0], box[0, 0, 2],
                   box[0, 1, 2], box[0, 2, 0]]  # convert to .gro format

    with open(args.top, 'r') as f:
        top = []
        for line in f:
            top.append(line)

    solvate_tails.add_water_placeholder(top)  # add a placeholder field to topol.top

    # choose whether we will place a water in the pores or in the tails
    placement_region = np.random.choice(['tails', 'pores'], p=[1 - args.probability, args.probability])

    if placement_region == 'pores':

        # find the pore centers
        pore_atom_indices = [a.index for a in t.topology.atoms if a.name in args.pore_atoms]
        pore_centers = physical.avg_pore_loc(4, t.xyz[0, pore_atom_indices, :])

        # randomly choose a pore from a uniform distribution
        pore = np.random.choice([0, 1, 2, 3])
        random_cylinder_point = random_pt_cylinder(args.radius_pore, t.unitcell_vectors[0, 2, 2])  # height is z-box vector
        water_placement_point = random_cylinder_point + [pore_centers[0, pore], pore_centers[1, pore], 0]  # shift so cylinder is centered at chosen pore center

    else:

        tail_atom_indices = [a.index for a in t.topology.atoms if a.name in args.tail_atoms]  # indices of oxygens in tails
        placement_atom = np.random.choice(tail_atom_indices)  # randomly choose which atom to place water molecule near
        # place point within spherical shell centered at placement_atom with inner radius=rmin and outside radius=rmax
        water_placement_point = solvate_tails.random_pt_spherical_shell(t.xyz[0, placement_atom, :], args.rmin, args.rmax)

    nrg = place_water(t.xyz[0, :, :], water_placement_point, ids, res, box_gromacs, args.emsteps, args.rem)

    while nrg > -50000:  # make sure the system is somewhat stable by checking potential energy after minimization
        nrg = place_water(t.xyz[0, :, :], water_placement_point, ids, res, box_gromacs, args.emsteps, args.rem)