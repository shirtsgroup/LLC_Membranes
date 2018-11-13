#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.path as path
import math
from llclib import physical
from llclib import transform
from llclib import file_rw
import os
from scipy import spatial
import place_solutes


def initialize():

    parser = argparse.ArgumentParser(description='Place solutes in the pores')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file where solutes will be placed')
    parser.add_argument('-r', '--ref', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Reference atom (s)'
                        ' used to locate pores')
    parser.add_argument('-n', '--number', default=1, type=int, help='Number of solutes to place in each pore')
    parser.add_argument('-p', '--top', help='Name of topology to modify. If not specified, dont update topology')
    parser.add_argument('-s', '--solute', default='ETH', help='Name of solute coordinate file to '
                        'place')
    parser.add_argument('--solute_path', default='/home/bcoscia/PycharmProjects/GitHub/HII/top/topologies', help='path '
                        'to where solute coordinate files are held')
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of layers')
    parser.add_argument('-nn', '--neighbors', default=3, type=int, help='Number of nearest neighbor water molecules to remove')
    parser.add_argument('-o', '--output', default='solutes.gro', type=str, help='Name of output .gro file')

    args = parser.parse_args()

    return args


def put_in_box(pt, x_box, y_box, m, angle):
    """
    :param pt: The point to place back in the box
    :param x_box: length of box in x dimension
    :param y_box: length of box in y dimension
    :param m: slope of box vector
    :param angle: angle between x axis and y box vector
    :return: coordinate shifted into box
    """

    b = - m * x_box  # y intercept of box vector that does not pass through origin (right side of box)
    if pt[1] > m*pt[0]:  # if the point is on the left side of the box
        pt[0] += x_box
    if pt[1] < m*(pt[0] - b):  # if the point is on the right side of the box
        pt[0] -= x_box
    if pt[1] < 0:
        pt[:2] += [np.cos(angle)*x_box, np.sin(angle)*x_box]  # if the point is under the box
    if pt[1] > y_box:
        pt[:2] -= [np.cos(angle)*x_box, np.sin(angle)*x_box]

    return pt


def trace_pores(pos, box, layers):
    """
    Find the line which traces through the center of the pores
    :param pos: positions of atoms used to define pore location (args.ref) [natoms, 3]
    :param box: xy box vectors, [2, 2], mdtraj format
    :param layers: number of layers
    :return: points which trace the pore center
    """

    atoms_p_pore = int(pos.shape[0] / 4)  # atoms in each pore
    atoms_p_layer = int(atoms_p_pore / layers)  # atom per layer

    v = np.zeros([4, 2])  # vertices of unitcell box
    v[0, :] = [0, 0]
    v[1, :] = [box[0, 0], 0]
    v[3, :] = [box[1, 0], box[1, 1]]
    v[2, :] = v[3, :] + [box[0, 0], 0]

    center = [np.mean(v[:, 0]), np.mean(v[:, 1]), 0]  # geometric center of box
    bounds = path.Path(v)  # create a path tracing the vertices, v

    angle = np.arccos(box[1, 1]/box[0, 0])  # angle of monoclinic box
    if box[1, 0] < 0:  # the case of an obtuse angle
        angle += np.pi / 2

    m = (v[3, 1] - v[0, 1]) / (v[3, 0] - v[0, 0])  # slope from points connecting first and third vertices

    centers = np.zeros([4*layers, 3])

    for p in range(4):
        pore = pos[p*atoms_p_pore:(p+1)*atoms_p_pore, :]  # coordinates for atoms belonging to a single pore
        for l in range(layers):
            before = pore[l*atoms_p_layer, :]  # choose the first atom as a reference

            shift = transform.translate(pore[l*atoms_p_layer:(l+1)*atoms_p_layer, :], before, center)  # shift everything to towards the center

            for i in range(shift.shape[0]):  # check if the points are within the bounds of the unitcell
                if not bounds.contains_point(shift[i, :2]):
                    shift[i, :] = put_in_box(shift[i, :], box[0, 0], box[1, 1], m, angle)  # if its not in the unitcell, shift it so it is

            c = np.zeros([1, 3])
            c[0, :] = [np.mean(shift[:, 0]), np.mean(shift[:, 1]), np.mean(shift[:, 2])]  # geometric center of reference atoms in this layer

            centers[p*layers + l, :] = transform.translate(c, center, before)  # move everything back to where it was

            if not bounds.contains_point(centers[p*layers, :]):  # make sure everything is in the box again
                centers[p*layers + l, :] = put_in_box(centers[p*layers + l, :], box[0, 0], box[1, 1], m, angle)

        #     print(centers[p*layers + l, :])
        # exit()

    return centers


def placement(z, pts, box):
    """
    :param z: z location where solute should be placed
    :param pts: points which run through the pore
    :return: location to place solute
    """

    v = np.zeros([4, 2])  # vertices of unitcell box
    v[0, :] = [0, 0]
    v[1, :] = [box[0, 0], 0]
    v[3, :] = [box[1, 0], box[1, 1]]
    v[2, :] = v[3, :] + [box[0, 0], 0]
    center = [np.mean(v[:, 0]), np.mean(v[:, 1]), 0]  # geometric center of box
    bounds = path.Path(v)  # create a path tracing the vertices, v

    angle = np.arccos(box[1, 1]/box[0, 0])  # angle of monoclinic box
    if box[1, 0] < 0:  # the case of an obtuse angle
        angle += np.pi / 2

    m = (v[3, 1] - v[0, 1]) / (v[3, 0] - v[0, 0])  # slope from points connecting first and fourth vertices

    # shift = transform.translate(z, before, center)
    #
    # put_in_box(pt, box[0, 0], box[1, 1], m, angle)

    # find z positions, in between which solute will be placed
    lower = 0
    while pts[lower, 2] < z:
        lower += 1

    upper = pts.shape[0] - 1
    while pts[upper, 2] > z:
        upper -= 1

    limits = np.zeros([2, 3])
    limits[0, :] = pts[lower, :]
    limits[1, :] = pts[upper, :]

    shift = transform.translate(limits, limits[0, :], center)  # shift limits to geometric center of unit cell

    for i in range(shift.shape[0]):  # check if the points are within the bounds of the unitcell
        if not bounds.contains_point(shift[i, :2]):
            shift[i, :] = put_in_box(shift[i, :], box[0, 0], box[1, 1], m, angle)

    # Use parametric representation of line between upper and lower points to find the xy value where z is satsified
    v = shift[1, :] - shift[0, :]  # direction vector

    t = (0 - shift[0, 2]) / v[2]  # solve for t since we know z
    x = shift[0, 0] + t*v[0]
    y = shift[0, 1] + t*v[1]

    place = np.zeros([1, 3])
    place[0, :] = [x, y, 0]
    place = transform.translate(place, center, limits[0, :])

    if not bounds.contains_point(place[0, :]):  # make sure everything is in the box again
        place[0, :] = put_in_box(place[0, :], box[0, 0], box[1, 1], m, angle)

    return place


# class System(object):
#
#     def __init__(self, initial, solute, reference_atoms, layers):
#
#         t = md.load(initial)
#         self.xyz = t.xyz[0, :, :]
#         self.full_box = t.unitcell_vectors
#         self.box = t.unitcell_vectors[0, :2, :2]  # get xy box vectors
#         self.zmax, self.zmin = physical.thickness(initial, reference_atoms, grid=False)[1:3]  # limits for solute placement
#         self.layers = layers
#         self.angle = np.arccos(self.box[1, 1] / self.box[0, 0])
#         if self.box[1, 0] < 0:
#             self.angle += np.pi / 2
#
#         # get indices of reference atoms which will be used to calculate the location of the pores
#         keep = [a.index for a in t.topology.atoms if a.name in reference_atoms]
#         self.reference_positions = self.xyz[keep, :]  # get positions of reference atoms
#
#         # create spline that traces pore center
#         self.centers = trace_pores(self.reference_positions, self.box, layers)
#
#         # create an object for the solute to be placed
#         self.solute = place_solutes.Solute(solute)
#
#     def add_solutes(self, n):
#
#         """
#         :param n: number of solutes to add
#         """
#
#         pores = int(self.centers.shape[0] / self.layers)
#         solute_locations = np.zeros([pores * n, 3])
#         solute_coords = np.zeros([self.solute.xyz.shape[1] * n * pores, 3])
#         for i in range(pores):
#             zmin = np.min(self.centers[i * self.layers:(i + 1) * self.layers, 2])
#             zmax = np.max(self.centers[i * self.layers:(i + 1) * self.layers, 2])
#             z_placement = np.linspace(zmin, zmax, n + 2)[1:-1]  # points along z axis where solutes will be placed (exclude end points ... for now)
#             for j in range(args.number):
#                 place = placement(z_placement[j], self.centers[i * self.layers:(i + 1) * self.layers, :], self.box)
#                 solute_locations[i * n + j, :] = place
#                 solute_coords[
#                 (i * args.number + j) * self.solute.xyz.shape[1]:(i * n + j + 1) * self.solute.xyz.shape[1], :] = \
#                     transform.translate(self.solute.xyz[0, :, :], self.com, place[0, :])


if __name__ == "__main__":

    args = initialize()

    # TODO in the future : class-ify
    # sys = System(args.gro, args.solute, args.ref, args.layers)
    # sys.add_solutes(args.number)

    t = md.load(args.gro)  # load trajectory
    xyz = t.xyz[0, :, :]  # positions of all atoms

    full_box = t.unitcell_vectors
    box = t.unitcell_vectors[0, :2, :2]  # get xy box vectors
    zmax, zmin = physical.thickness(args.gro, args.ref, grid=False)[1:3]  # find limits for solute placement

    angle = np.arccos(box[1, 1]/box[0, 0])
    if box[1, 0] < 0:
        angle += np.pi / 2

    keep = [a.index for a in t.topology.atoms if a.name in args.ref]  # keep reference atoms

    pos = t.atom_slice(keep).xyz[0, :, :]  # get positions of reference atoms

    #Now find the xy locations where solutes should be placed (i.e. find the pores)
    centers = trace_pores(pos, box, args.layers)

    #Place solutes molecules by their center of mass
    solute = md.load('%s/%s.gro' % (args.solute_path, args.solute))

    solute_atom_names = [a.name for a in solute.topology.atoms]
    solute_resnames = [a.residue.name for a in solute.topology.atoms]

    com = np.zeros([3])  # center of mass of solute
    for i in range(solute.xyz.shape[1]):
        com += solute.xyz[0, i, :]
    com /= solute.xyz.shape[1]

    pores = int(centers.shape[0] / args.layers)
    solute_locations = np.zeros([pores*args.number, 3])
    solute_coords = np.zeros([solute.xyz.shape[1]*args.number*pores, 3])
    for i in range(pores):
        zmin = np.min(centers[i*args.layers:(i+1)*args.layers, 2])
        zmax = np.max(centers[i*args.layers:(i+1)*args.layers, 2])
        z_placement = np.linspace(zmin, zmax, args.number + 2)[1:-1]  # points along z axis where solutes will be placed (exclude end points ... for now)
        for j in range(args.number):
            place = placement(z_placement[j], centers[i*args.layers:(i+1)*args.layers, :], box)
            solute_locations[i*args.number + j, :] = place
            solute_coords[(i*args.number + j)*solute.xyz.shape[1]:(i*args.number + j + 1)*solute.xyz.shape[1], :] = \
                transform.translate(solute.xyz[0, :, :], com, place[0, :])
            # xyz = np.concatenate((xyz, transform.translate(solute.xyz[0, :, :], com, place[0, :])))

    # remove waters close to inserted solutes
    water = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']

    tree = spatial.cKDTree(xyz[water, :])
    rm = []

    for i in range(solute_locations.shape[0]):
        for j in tree.query(solute_locations[i, :], k=args.neighbors)[1]:
            rm.append(water[j])
            rm.append(water[j] + 1)
            rm.append(water[j] + 2)

    keep = [a.index for a in t.topology.atoms if a.residue.name != 'HOH']
    for i in water:
        if i not in rm:
            keep.append(i)
            keep.append(i + 1)
            keep.append(i + 2)

    pos = np.concatenate((xyz[keep, :], solute_coords))

    names = [a.name for a in t.topology.atoms if a.index not in rm]  # names of all atoms
    res = [a.residue.name for a in t.topology.atoms if a.index not in rm]  # residue names of all atoms

    # change names and residues back to be consistent with tip3p water model. idk why mdtraj does this to me D':
    for i, name in enumerate(res):
        if name == "HOH":
            res[i] = 'SOL'
            if names[i] == 'O':
                names[i] = 'OW'
            elif names[i] == 'H1':
                names[i] = 'HW1'
            elif names[i] == 'H2':
                names[i] = 'HW2'

    ids = names + solute_atom_names*args.number*4
    res_names = res + solute_resnames*args.number*4

    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                   full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]

    file_rw.write_gro_pos(pos, '%s' % args.output, box=box_gromacs, ids=ids, res=res_names)

    # update topology
    if args.top:

        topology = []
        with open('%s' % args.top, 'r') as f:
            for line in f:
                topology.append(line)

        system_index = 0
        while topology[system_index].count('[ system ]') == 0:
            system_index += 1

        topology.insert(system_index, ';%s Topology\n' % args.solute)
        topology.insert(system_index + 1, '#include "%s/%s.itp"\n' % (args.solute_path, args.solute))
        topology.insert(system_index + 2, '\n')

        # modify number of water molecules since some were removed
        sol = 0
        while topology[sol].count('SOL') == 0:
            sol += 1

        topology[sol] = '%s                 %d\n' % ('SOL', int(int(topology[sol].split()[1]) - args.neighbors*args.number*pores))
        topology.append('%s                 %s' % (solute_resnames[0], args.number*4))

        with open('solute.top', 'w') as f:
            for line in topology:
                f.write(line)