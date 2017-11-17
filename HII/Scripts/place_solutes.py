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


def initialize():

    parser = argparse.ArgumentParser(description='Place solutes in the pores')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file where solutes will be placed')
    parser.add_argument('-r', '--ref', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Reference atom (s)'
                        ' used to locate pores')
    parser.add_argument('-n', '--number', default=1, type=int, help='Number of solutes to place in each pore')
    parser.add_argument('-p', '--top', default='topol.top', help='Name of topology to modify')
    parser.add_argument('--ngrid', default=3, type=float, help='Number of grid points in one dimension')
    parser.add_argument('-s', '--solute', default='ethanol', help='Name of solute coordinate file to '
                        'place')
    parser.add_argument('--solute_path', default='/home/bcoscia/PycharmProjects/GitHub/HII/top/solutes', help='path to '
                        'where solute coordinate files are held')
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of layers')

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

    if pt[1] < m*pt[0]:  # if the point is on the left side of the box
        pt[0] += x_box
    if pt[1] > m*(pt[0] - x_box):  # if the point is on the right side of the box
        pt[0] -= x_box
    if pt[1] < 0:
        pt[:2] += [np.cos(angle)*box[0, 0], np.sin(angle)*box[0, 0]]  # if the point is under the box
    if pt[1] > y_box:
        pt[:2] -= [np.cos(angle)*box[0, 0], np.sin(angle)*box[0, 0]]

    return pt


def trace_pores(pos, box, layers):
    """
    Find the line which traces through the center of the pores
    :param pos: positions of atoms used to define pore location (args.ref) [natoms, 3]
    :param box: xy box vectors, [2, 2], mdtraj format
    :return: points which trace the pore center
    """

    atoms_p_pore = int(pos.shape[0] / 4)
    atoms_p_layer = int(atoms_p_pore / layers)

    v = np.zeros([4, 2])  # vertices of unitcell box
    v[0, :] = [0, 0]
    v[1, :] = [box[0, 0], 0]
    v[3, :] = [box[1, 0], box[1, 1]]
    v[2, :] = v[3, :] + [box[0, 0], 0]
    center = [np.mean(v[:, 0]), np.mean(v[:, 1]), 0]
    bounds = path.Path(v)

    angle = np.arccos(box[1, 1]/box[0, 0])
    if box[1, 0] < 0:
        angle += np.pi / 2

    m = (v[3, 1] - v[0, 1]) / (v[3, 0] - v[0, 0])  # slope from points connecting first and third vertices

    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                   full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]

    centers = np.zeros([4*layers, 3])

    for p in range(4):
        pore = pos[p*atoms_p_pore:(p+1)*atoms_p_pore, :]  # coordinates for atoms belonging to a single pore
        for l in range(layers):
            before = pore[l*atoms_p_layer, :]  # choose the first atom as a reference
            shift = transform.translate(pore[l*atoms_p_layer:(l+1)*atoms_p_layer, :], before, center)  # shift everything to towards the center

            for i in range(shift.shape[0]):
                if not bounds.contains_point(shift[i, :2]):
                    shift[i, :] = put_in_box(shift[i, :], box[0, 0], box[1, 1], m, angle)

            c = np.zeros([1, 3])
            c[0, :] = [np.mean(shift[:, 0]), np.mean(shift[:, 1]), np.mean(shift[:, 2])]
            centers[p*layers + l, :] = transform.translate(c, center, before)

            if not bounds.contains_point(centers[p*layers, :]):
                centers[p*layers + l, :] = put_in_box(centers[p*layers, :], box[0, 0], box[1, 1], m, angle)

    return centers


def placement(z, pts):
    """
    :param z: z location where solute should be placed
    :param pts: points which run through the pore
    :return: location to place solute
    """

    lower = 0
    while pts[lower, 2] < z:
        lower += 1

    upper = pts.shape[0] - 1
    while pts[upper, 2] > z:
        upper -= 1

    # Use parametric representation of line between upper and lower points to find the xy value where z is satsified

    v = pts[upper, :] - pts[lower, :]  # direction vector
    t = (z - pts[lower, 2]) / v[2]  # solve for t since we know z
    x = pts[lower, 0] + t*v[0]
    y = pts[lower, 1] + t*v[1]

    return [x, y, z]


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)  # load trajectory
    xyz = t.xyz[0, :, :]

    names = [a.name for a in t.topology.atoms]
    res = [a.residue.name for a in t.topology.atoms]

    # change names and residues back to be consistent with tip3p water model. idk mdtraj does this to me D':
    for i, name in enumerate(res):
        if name == "HOH":
            res[i] = 'SOL'
            if names[i] == 'O':
                names[i] = 'OW'
            elif names[i] == 'H1':
                names[i] = 'HW1'
            elif names[i] == 'H2':
                names[i] = 'HW2'

    full_box = t.unitcell_vectors
    box = t.unitcell_vectors[0, :2, :2]  # get xy box vectors
    zmax, zmin = physical.thickness(args.gro, args.ref, grid=False)[1:3]  # find limits for solute placement
    z_placement = np.linspace(zmin, zmax, args.number + 2)[1:-1]  # points along z axis where solutes will be placed (exclude end points ... for now)

    angle = np.arccos(box[1, 1]/box[0, 0])
    if box[1, 0] < 0:
        angle += np.pi / 2

    keep = [a.index for a in t.topology.atoms if a.name in args.ref]  # keep reference atoms

    pos = t.atom_slice(keep).xyz[0, :, :]  # get positions of reference atoms

    # Now find the xy locations where solutes should be placed (i.e. find the pores)
    centers = trace_pores(pos, box, args.layers)

    # Place solutes molecules by their center of mass
    solute = md.load('%s/%s.gro' % (args.solute_path, args.solute))

    solute_atom_names = [a.name for a in solute.topology.atoms]
    solute_resnames = [a.residue.name for a in solute.topology.atoms]

    com = np.zeros([3])  # center of mass of solute
    for i in range(solute.xyz.shape[1]):
        com += solute.xyz[0, i, :]

    com /= solute.xyz.shape[1]
    pores = int(centers.shape[0] / args.layers)
    for i in range(pores):
        for j in range(args.number):
            place = placement(z_placement[j], centers[i*args.layers:(i+1)*args.layers, :])
            xyz = np.concatenate((xyz, transform.translate(solute.xyz[0, :, :], com, place)))

    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                   full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]

    ids = names + solute_atom_names*args.number*4
    res_names = res + solute_resnames*args.number*4

    file_rw.write_gro_pos(xyz, 'solutes.gro', box=box_gromacs, ids=ids, res=res_names)

    # update topology

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
    topology.append('%s                 %s' % (solute_resnames[0], args.number*4))

    with open('solute_top.top', 'w') as f:
        for line in topology:
            f.write(line)