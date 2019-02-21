#!/usr/bin/env python
import numpy as np
import math


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


def rotate_vector(xyz, v1, v2):
    """
    :param xyz: xyz coordinates of object to be rotated
    :param v1: original vector
    :param v2: direction you want v1 to be pointing in
    :return: rotated coordinates
    """

    quad = quadrant(v1)
    # first find the angle between v1 and v2
    num = np.dot(v1, v2)
    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    theta = np.arccos(num / denom)

    if quad == 1 or quad == 2:
        Rz = rotate_z(-theta)
    else:
        Rz = rotate_z(theta)

    pos = np.zeros_like(xyz)
    for i in range(np.shape(xyz)[0]):
        pos[i, :] = np.dot(Rz, xyz[i, :])

    return pos