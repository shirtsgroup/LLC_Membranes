#!/usr/bin/python
"""
Reposition monomers to surround the average pore center as a way of aiding equilibration
"""

import argparse
import Structure_char
import restrain
import math
import numpy as np
import build


def initialize():
    parser = argparse.ArgumentParser(description='Build LLC Structure')

    parser.add_argument('-f', '--file', default='wiggle.gro', type=str, help='.gro coordinate input file')
    parser.add_argument('-r', '--radius', default=.2, type=float, help='distance from pore center to reposition monomers')
    parser.add_argument('-a', '--atoms', default=['C6', 'C'], help='Atoms which will be used to determine pore centers.'
                                                                   'The first entry in the list should be the atom '
                                                                   'which you want to rotate the monomer about')
    args = parser.parse_args()

    return args


def slope(pt1, pt2):
    """
    Given two points, find the 2D slope
    :param pt1: xy coordinates of 1st point
    :param pt2: xy coordinates of 2nd point
    :return: slope between points
    """
    m = (pt1[1] - pt2[1])/(pt1[0] - pt2[0])  # slope
    m = float(m)
    return m


def angles(pts, pcenters):
    """
    Given 2 points, find the angle needed to rotate the monomer so that it is facing the center of the pore. NOTE: this
    will only use the xy coordinates and rotate with respect to the xy plane
    :param pts: 3 xyz coordinates which should roughly form a line
    :param pcenters: the xy locations of the pore centers
    :return: the angle needed to rotate each monomer
    """
    nMon = pts.shape[1] / 2  # divide by 2 since there are two points for each monomer
    npores = pcenters.shape[1]
    monppore = nMon / npores
    thetas = np.zeros([nMon])

    for i in range(npores):
        for j in range(monppore):
            pt1 = pts[:2, i*monppore + 2*j]  # first entry in args.atoms
            pt2 = pts[:2, i*monppore + 2*j + 1]  # second entry in args.atoms
            m1 = slope(pt1, pt2)
            pt3 = pcenters[:2, i]
            m2 = slope(pt1, pt3)
            thetas[i*monppore + j] = -math.atan((m1 - m2)/(1 + m1*m2))  # find angle between lines

    return thetas


def translate(obj, pt=np.array([0, 0, 0]), r='na'):
    """
    Create a translation matrix to translate an object by a vector v
    :param obj: the position of the object being translated
    :param pt: the point to which you want to translate the object (optional, default is the origin)
    :param r: distance to translate an object at fixed angle. Useful for translating a point radially outward from the
    origin
    :return: translation matrix
    """
    if r == 'na':
        T = np.zeros([4, 4])  # Create the skeleton of a translation matrix
        for i in range(4):
            for j in range(4):
                if i == j:
                    T[i, j] = 1  # the diagonal is all one

        v = pt - obj  # direction in which to move the object

        T[:2, 3] = v[:2]  # the last column contains the vector in whose direction we wish to translate the object. We
                          # also do not want to translate in the z direction in this case

    else:
        vx, vy = build.transdir(obj, origin=pt)
        theta = math.atan(obj[1]/obj[0])

        T = np.zeros([4, 4])  # Create the skeleton of a translation matrix
        for i in range(4):
            for j in range(4):
                if i == j:
                    T[i, j] = 1  # the diagonal is all one

        T[0, 3] = vx*r*math.cos(theta)
        T[1, 3] = vy*r*math.sin(theta)

    return T


def quadrant(pt):  # looks at [x,y] values and determines which quadrant the point is in
    if pt[0] > 0 and pt[1] > 0:
        return 1
    elif pt[0] < 0 and pt[1] < 0:
        return 3
    elif pt[1] < 0 < pt[0]:
        return 4
    elif pt[0] < 0 < pt[1]:
        return 2
    else:
        return 0  # the case where the point lies on the x or y axis


def rmatrix(theta):
    """
    Create a rotation matrix to rotate points w.r.t. the xy plane given an angle, theta
    :param theta: angle of rotation - radians
    :param pt: the point around which to rotate
    :return: rotation matrix
    """
    Rx = np.zeros([3, 3])  # makes a 3 x 3 zero matrix
    Rx[0, 0] = math.cos(theta)  # This line and subsequent edits to Rx fills in entries needed for rotation matrix
    Rx[1, 0] = math.sin(theta)
    Rx[0, 1] = -math.sin(theta)
    Rx[1, 1] = math.cos(theta)
    Rx[2, 2] = 1

    return Rx


def reposition(pts, nMon, pcenters, atom_no, t_atom, valence):
    """
    Rotate monomers to face toward the pore centers
    :param pts: all of the coordinates from coordinate file in a numpy array: [xyz, n_atoms]
    :param thetas: angles by which to rotate each monomer
    :return: all of the monomers rotated in the correct direction
    """
    atoms = pts.shape[1]
    atomspmon = atoms / nMon - 1  # atoms per monomer
    t_origin = np.zeros(pts.shape)  # translate to origin
    rotated = np.zeros(pts.shape)  # rotate
    t_back = np.zeros(pts.shape)  # translate back to starting point
    pores = pcenters.shape[1]
    monppore = pts.shape[1] / pores

    # for i in range(nMon):
    #     loc = pts[:, i*atomspmon + atom_no]  # the xy coordinates of the chosen atom of this monomer
    #     T_fwd = translate(loc)  # translate from loc to origin
    #     T_bwd = translate(-loc)  # translate from origin back to loc
    #     Rx = rmatrix(thetas[i])  # rotation matrix at origin
    #     for j in range(atomspmon):
    #         cat = np.append(pts[:, i*atomspmon + j], [1])
    #         t_origin[:, i*atomspmon + j] = np.dot(T_fwd, cat)[:3]
    #     for j in range(atomspmon):
    #         cat = t_origin[:, i*atomspmon + j]
    #         rotated[:, i*atomspmon + j] = np.dot(Rx, cat)[:3]
    #     for j in range(atomspmon):
    #         cat = np.append(rotated[:, i*atomspmon + j], [1])
    #         t_back[:, i*atomspmon + j] = np.dot(T_bwd, cat)[:3]

    for i in range(nMon):
        pore = int(i / (nMon/pores))
        center = np.append(pcenters[:, pore], [0])
        loc = pts[:, i*atomspmon + atom_no]
        T_fwd = translate(loc, pt=center)  # translate from loc to origin
        for j in range(atomspmon):
            cat = np.append(pts[:, i*atomspmon + j], [1])
            t_origin[:, i*atomspmon + j] = np.dot(T_fwd, cat)[:3]
        cat = np.append(pts[:, atoms - nMon + i], [1])
        t_origin[:, atoms - nMon + i] = np.dot(T_fwd, cat)[:3]
        T_bwd = translate(t_origin[:, i*atomspmon + t_atom], r=args.radius, pt=center)
        for j in range(atomspmon):
            cat = np.append(t_origin[:, i*atomspmon + j], [1])
            t_back[:, i*atomspmon + j] = np.dot(T_bwd, cat)[:3]
        cat = np.append(t_origin[:, atoms - nMon + i], [1])
        t_back[:, atoms - nMon + i] = np.dot(T_bwd, cat)[:3]

    return t_back


def sort(pts, nMon, valence):

    atoms = pts.shape[1]
    atomspmon = pts.shape[1] / nMon

    sorted_coords = np.zeros(pts.shape)
    for i in range(nMon):
        for j in range(atomspmon - valence):
            sorted_coords[:, i*(atomspmon - valence) + j] = pts[:, i*atomspmon + j]
        sorted_coords[:, atoms - nMon*valence + i] = pts[:, i*atomspmon + j + 1]

    return sorted_coords


def write_gro(pts, ids, res, name="repositioned.gro"):
    """
    Write a .gro file given coordinates, atomic identities and which residue they belong to
    :param pts: atomic coordinates, np array [xyz, n_atoms]
    :param ids: identities, np array [n_atoms]
    :param res: residue names, np array [n_atoms]
    :return:
    """

    atoms = ids.shape[0]
    f = open(name, 'w')
    f.write('This is a .gro file\n')
    f.write('%s\n' % atoms)

    n_res = 1
    for i in range(atoms):
        f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(n_res, '%s' % res[i], '%s' % ids[i], i + 1, pts[0, i],
                                                                       pts[1, i], pts[2, i]))
        if i != atoms - 1 and res[i + 1] != res[i] or res[i] == 'NA':
            n_res += 1

    f.write('0.0000000 0.0000000 0.0000000\n')
    f.close()

if __name__ == "__main__":

    args = initialize()

    f = open(args.file, 'r')
    gro = []
    for line in f:
        gro.append(line)
    f.close()

    coords, atoms, all_coords, ids, res = restrain.get_coordinates(gro, 2, 'NA')

    atom_no = np.where(ids == args.atoms[0])[0][0]  # meh
    t_atom = np.where(ids == args.atoms[1])[0][0]  # double meh

    p_centers = Structure_char.avg_pore_loc(4, coords, atoms)

    thetas = angles(coords, p_centers)

    translated = reposition(all_coords, coords.shape[1], p_centers, atom_no, t_atom, 1)

    # repositioned = sort(translated, coords.shape[1], 1)

    write_gro(translated, ids, res)