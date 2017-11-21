#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import lc_class
import matplotlib.pyplot as plt
import matplotlib.path as path
import math
from llclib import physical
from llclib import transform
from llclib import file_rw
import os
from scipy import spatial


def initialize():

    parser = argparse.ArgumentParser(description='Straighten the pores of a solvates system')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file where solutes will be placed')
    parser.add_argument('-b', '--build_monomer', default='NAcarb11Vd.gro', help='Name of monomer used to build system')
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
        pt[:2] += [np.cos(angle)*x_box, np.sin(angle)*x_box]  # if the point is under the box
    if pt[1] > y_box:
        pt[:2] -= [np.cos(angle)*x_box, np.sin(angle)*x_box]

    return pt


def shift_pores(pos, box, layers):
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

if __name__ == "__main__":

    args = initialize()

    t = md.load('%s' % args.gro)
    mon = lc_class.LC('%s' % args.build_monomer)

    monomer = [a.index for a in t.topology.atoms if a.residue.name != 'HOH']
    nmon = len(monomer) // (mon.natoms + mon.no_ions)
    mon_per_pore = nmon / 4