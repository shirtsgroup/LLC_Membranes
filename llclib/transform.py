#! /usr/bin/env python

"""
Coordinate transforms and related manipulations of positions
"""

import numpy as np
import math


def layer_dist(layers, nopores, distribution, monomers, alt_1, alt_2):

    layer_distribution = np.zeros([layers*nopores], dtype=int)

    if distribution == 'uniform':
        for i in range(0, len(layer_distribution)):
            layer_distribution[i] = monomers
    if distribution == 'alternating':
        for i in range(layer_distribution.shape[-1]):
            if i % 2 == 0:
                layer_distribution[i] = alt_1
            if i % 2 == 1:
                layer_distribution[i] = alt_2
    return layer_distribution


def slope(pt1, pt2):
    m = (pt1[1] - pt2[1])/(pt1[0] - pt2[0])  # slope
    return m


def rotateplane(plane, angle=0):
    """
    Calculate a rotation matrix to rotate a plane in 3 dimensions
    :param plane: indices of atoms making up plane which is being aligned in the xy plane
    :param angle: desired angle between xy plane (optional, default = 0 i.e. in plane)
    :return:
    """

    # vector pointing from point 1 to point 2
    v12 = plane[1, :] - plane[0, :]
    v13 = plane[2, :] - plane[0, :]

    # The cross product of v12 and v13 give a vector that is perpendicular to the plane:
    N = np.cross(v12, v13)
    N_desired = [0, math.sin(angle), math.cos(angle)]  # vector in the direction normal to our desired plane orientation

    RotationAxis = np.cross(N, N_desired)

    theta = math.acos(np.dot(N, N_desired)/(np.linalg.norm(N)*np.linalg.norm(N_desired)))  #  Rotation Angle (radians)

    L = [RotationAxis[0]/np.linalg.norm(RotationAxis), RotationAxis[1]/np.linalg.norm(RotationAxis),
                       RotationAxis[2]/np.linalg.norm(RotationAxis)]  # normalized Rotation Axis
    # ^ see: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/

    u, v, w = L[0], L[1], L[2]

    R = np.zeros((4, 4))
    R[3, 3] = 1
    R[0, 0] = u**2 + (v**2 + w**2)*math.cos(theta)  # math.cos takes theta in radians by default
    R[0, 1] = u*v*(1 - math.cos(theta)) - w*math.sin(theta)
    R[0, 2] = u*w*(1 - math.cos(theta)) + v*math.sin(theta)
    R[1, 0] = u*v*(1 - math.cos(theta)) + w*math.sin(theta)
    R[1, 1] = v**2 + (u**2 + w**2)*math.cos(theta)
    R[1, 2] = v*w*(1 - math.cos(theta)) - u*math.sin(theta)
    R[2, 0] = u*w*(1 - math.cos(theta)) - v*math.sin(theta)  # math.cos takes theta in radians by default
    R[2, 1] = w*v*(1 - math.cos(theta)) + u*math.sin(theta)
    R[2, 2] = w**2 + (u**2 + v**2)*math.cos(theta)

    return R


def quadrant(pt, origin=[0, 0]):
    """ Find out which quadrant of the xy plane a point is sitting in
    II    |    I
          |
    -------------
          |
    III   |    IV
    :param: pt: point to be tested
    :param: origin: the location of the origin. Default is [0, 0] but can be set arbitrarily (such as a pore center)
    """
    if pt[0] > origin[0] and pt[1] > origin[1]:
        return 1
    elif pt[0] < origin[0] and pt[1] < origin[1]:
        return 3
    elif pt[0] > origin[0] and pt[1] < origin[1]:
        return 4
    elif pt[0] < origin[0] and pt[1] > origin[1]:
        return 2
    else:
        return 0  # the case where the point lies on the x or y axis


def transdir(pt, origin=[0, 0]):
    """
    figure out in which direction the coordinates will be shifted. They are always shifted away from the origin
    :param pt:
    :return:
    """
    if quadrant(pt, origin) == 1:  # e.g. in quadrant 1, the x's are shifted in the positive x and positive y directions
        vx = 1
        vy = 1
    elif quadrant(pt, origin) == 2:  # in quadrant 2, the x's are shifted negative and the y's are shifted positive
        vx = -1
        vy = 1
    elif quadrant(pt, origin) == 3:  # in quadrant 3, the x's and y's are shifted down
        vx = -1
        vy = -1
    elif quadrant(pt, origin) == 4:  # in quadrant 4, the x's are shifted positive and the y's are shifted negative
        vx = 1
        vy = -1
    # These next three conditionals are very unlikely but are included for completeness and to avoid future errors
    elif quadrant(pt, origin) == 0 and pt[0] == origin[0]:  # i.e., it lies on the y - axis
        if pt[1] > 0:  # the point is on the positive y-axis
            vx = 0  # no x-shift
            vy = 1  # shift in the positive y direction
        if pt[1] < 0:  # the point is on the negative y-axis
            vx = 0  # no x-shift
            vy = -1  # shift in the negative y direction
    elif quadrant(pt, origin) == 0 and pt[1] == origin[1]:  # i.e., it lies on the x - axis
        if pt[0] > 0:  # the point is on the positive x-axis
            vx = 1  # shift in the positive x direction
            vy = 0  # no y-shift
        if pt[0] < 0:  # the point is on the negative y-axis
            vx = -1  # shift in the negative x direction
            vy = 0  # no y-shift

    return vx, vy


def rotate_z(theta):
    """
    formerly 'rotate'
    :param: angle by which to rotate the monomer
    :return: Rotation matrix to rotate input vector about z-axis
    """
    Rx = np.zeros([3, 3])  # makes a 3 x 3 zero matrix
    Rx[0, 0] = math.cos(theta)
    Rx[1, 0] = math.sin(theta)
    Rx[0, 1] = -math.sin(theta)
    Rx[1, 1] = math.cos(theta)
    Rx[2, 2] = 1

    return Rx


def reposition(xyz, R, ref_index, lineatoms, pore_radius):

    b = np.ones([1])
    for i in range(np.shape(xyz)[1]):
        coord = np.concatenate((xyz[:, i], b))
        x = np.dot(R, coord)
        xyz[:, i] = x[:3]

    # Now translate the structure to the origin

    translation = np.matrix([[1, 0, 0, -xyz[0, ref_index]], [0, 1, 0, -xyz[1, ref_index]],
                         [0, 0, 1, -xyz[2, ref_index]], [0, 0, 0, 1]])

    b = np.ones([1])
    for i in range(np.shape(xyz)[1]):
        coord = np.concatenate((xyz[:, i], b))
        x = np.dot(translation, coord)
        xyz[:, i] = x[0, :3]

    # Now rotate the xy coordinates so that the molecule is pointing towards the origin

    pt1 = [xyz[0, lineatoms[0]], xyz[1, lineatoms[0]]]  # location of C
    pt2 = [xyz[0, lineatoms[1]], xyz[1, lineatoms[1]]]  # location of C3

    origin = [0, 0]

    # find slope between two points

    m1 = slope(pt1, pt2)

    m2 = 0  # slope of line y = 0

    # find angle between lines

    theta = -math.atan((m1 - m2)/(1 + m1*m2))

    vx, vy = transdir(pt1)

    # Translation matrix
    translation = np.matrix([[1, 0, 0, vx*pore_radius*math.cos(theta)], [0, 1, 0, vy*pore_radius*math.sin(theta)],\
                             [0, 0, 1, 0], [0, 0, 0, 1]])

    b = np.ones([1])
    for i in range(np.shape(xyz)[1]):
        coord = np.concatenate((xyz[:, i], b))
        x = np.dot(translation, coord)
        xyz[:, i] = x[0, :3]

    return xyz