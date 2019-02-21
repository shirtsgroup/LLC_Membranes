#! /usr/bin/env python

"""
Commonly used coordinate transforms and related manipulations of positions, including fourier transforms
"""
from __future__ import division
from __future__ import print_function

import numpy as np
import math
import tqdm


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
    m = (pt1[1] - pt2[1]) / (pt1[0] - pt2[0])  # slope
    return m


def plane_rotation_matrix(n, angle):
    """ Calculate the rotation matrix required to rotate a plane in 3 dimensions

    :param n: vector normal to plane to be rotated
    :param angle:

    :return:
    """


def rotateplane(plane, angle=0):
    """
    Calculate a rotation matrix to rotate a plane in 3 dimensions
    :param plane: indices of atoms making up plane which is being aligned in the xy plane
    :param angle: desired angle between xy plane (optional, default = 0 i.e. in plane)
    :return:
    """

    if plane[0, 2] == plane[1, 2] and plane[1, 2] == plane[2, 2]:

        print('Planes are already coplanar')

        return False

    else:
        # vector pointing from point 1 to point 2
        v12 = plane[1, :] - plane[0, :]
        v13 = plane[2, :] - plane[0, :]

        # The cross product of v12 and v13 give a vector that is perpendicular to the plane:
        N = np.cross(v12, v13)
        N_desired = [0, math.sin(angle), math.cos(angle)]  # vector in the direction normal to our desired plane orientation

        RotationAxis = np.cross(N, N_desired)

        theta = math.acos(np.dot(N, N_desired) / (np.linalg.norm(N)*np.linalg.norm(N_desired)))  #  Rotation Angle (radians)

        L = [RotationAxis[0] / np.linalg.norm(RotationAxis), RotationAxis[1] / np.linalg.norm(RotationAxis),
                           RotationAxis[2] / np.linalg.norm(RotationAxis)]  # normalized Rotation Axis
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


def rotateplane_coords(xyz, plane, angle=0):
    """
    Calculate a rotation matrix to rotate a plane in 3 dimensions
    :param xyz: xyz coordinates of all positions to be rotated
    :param plane: indices of atoms making up plane which is being aligned in the xy plane
    :param angle: desired angle between xy plane (optional, default = 0 i.e. in plane)
    :return:
    """

    if plane[0, 2] == plane[1, 2] and plane[1, 2] == plane[2, 2]:

        pass

    else:

        R = rotateplane(plane, angle=angle)

        b = np.ones([1])
        for i in range(np.shape(xyz)[0]):
            coord = np.concatenate((xyz[i, :], b))
            x = np.dot(R, coord)
            xyz[i, :] = x[:3]

    return xyz


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


def rotate_x(theta):
    """
    formerly 'rotate'
    :param: angle by which to rotate the monomer
    :return: Rotation matrix to rotate input vector about x-axis
    """
    Rx = np.zeros([3, 3])  # makes a 3 x 3 zero matrix
    Rx[0, 0] = 1
    Rx[2, 1] = math.sin(theta)
    Rx[1, 2] = -math.sin(theta)
    Rx[1, 1] = math.cos(theta)
    Rx[2, 2] = math.cos(theta)

    return Rx


def rotate_z(theta):
    """
    formerly 'rotate'
    :param: angle by which to rotate the monomer
    :return: Rotation matrix to rotate input vector about z-axis
    """
    Rz = np.zeros([3, 3])  # makes a 3 x 3 zero matrix
    Rz[0, 0] = math.cos(theta)
    Rz[1, 0] = math.sin(theta)
    Rz[0, 1] = -math.sin(theta)
    Rz[1, 1] = math.cos(theta)
    Rz[2, 2] = 1

    return Rz


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

    theta = -math.atan((m1 - m2) / (1 + m1*m2))

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


def translate(xyz, before, after):
    """
    :param xyz: coordinates of set of points to be translated [npts, 3]
    :param before: reference coordinate location before [3]
    :param after: reference coordinate location after [3]
    :return: translated points with respect to reference coordinate before/after locations [npts, 3]
    """

    pos = np.copy(xyz)
    direction = after - before

    translation = np.matrix([[1, 0, 0, direction[0]], [0, 1, 0, direction[1]],
                         [0, 0, 1, direction[2]], [0, 0, 0, 1]])

    b = np.ones([1])
    for i in range(pos.shape[0]):
        coord = np.concatenate((pos[i, :], b))
        x = np.dot(translation, coord)
        pos[i, :] = x[0, :3]

    return pos


def shift_matrices(images, angle, xbox, ybox):

    mat_dim = images * 2 + 1
    x_shift = np.zeros((mat_dim, mat_dim))
    y_shift = np.zeros((mat_dim, mat_dim))

    # shift in x and y direction due to angle
    x_comp = np.cos(angle*np.pi/180)*ybox
    y_comp = np.sin(angle*np.pi/180)*ybox

    for i in range(images + 1):
        x_shift[images, images + i] = i*xbox
        x_shift[images, images - i] = -i*xbox

    for i in range(1, images + 1):
        for j in range(images + 1):
            x_shift[images + i, images + j] = -i*x_comp + j*xbox
            x_shift[images - i, images + j] = i*x_comp + j*xbox
            x_shift[images + i, images - j] = -i*x_comp - j*xbox
            x_shift[images - i, images - j] = i*x_comp - j*xbox

    for i in range(1, images + 1):
        y_shift[images + i, :] = - i * y_comp
        y_shift[images - i, :] = i * y_comp

    return x_shift, y_shift


def pbcs(pts, images, angle, box, frame, nogap=False):
    """
    :param pts:
    :param images:
    :param angle:
    :param box:
    :param frame:
    :param nogap:
    :return:
    """

    xbox = np.linalg.norm(box[0, 0, :])
    ybox = np.linalg.norm(box[0, 1, :])
    zbox = np.linalg.norm(box[0, 2, :])

    x_shift, y_shift = shift_matrices(images, angle, xbox, ybox)
    mat_dim = 2 * images + 1

    tot_pts = np.shape(pts)[1]

    if nogap:
        translated_pts = np.zeros([3, 3*mat_dim**2, tot_pts])  # include all images in z direction as well
        z = [-1, 0, 1]
    else:
        translated_pts = np.zeros([3, mat_dim**2, tot_pts])
        z = [0]

    if len(pts.shape) == 3:
        for p in range(tot_pts):
            for k in range(len(z)):
                for i in range(mat_dim):
                    for j in range(mat_dim):
                        translated_pts[0, k*mat_dim**2 + i*mat_dim + j, p] = x_shift[i, j] + pts[frame, p, 0] # changed order for xlink.py and compatibility with mdtraj
                        translated_pts[1, k*mat_dim**2 + i*mat_dim + j, p] = y_shift[i, j] + pts[frame, p, 1]
                        translated_pts[2, k*mat_dim**2 + i*mat_dim + j, p] = pts[frame, p, 2] + zbox*z[k]
    else:
        for p in range(tot_pts):
            for i in range(mat_dim):
                for j in range(mat_dim):
                    translated_pts[0, i*mat_dim + j, p] = x_shift[i, j] + pts[p, 0] / 10
                    translated_pts[1, i*mat_dim + j, p] = y_shift[i, j] + pts[p, 1] / 10
                    translated_pts[2, i*mat_dim + j, p] = pts[p, 2] / 10  # z position unchanged

    return translated_pts


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


def rotate_coords_x(pos, angle):
    """
    :param xyz: xyz coordinates of atoms to be rotated, degrees
    :param angle: angle to rotate them by w.r.t origin
    :return:
    """

    xyz = np.copy(pos)
    angle *= (np.pi / 180)  # convert to radians

    R = rotate_x(angle)

    for i in range(np.shape(xyz)[0]):
        xyz[i, :] = np.dot(R, xyz[i, :])

    return xyz


def rotate_coords_z(pos, angle):
    """
    :param xyz: xyz coordinates of atoms to be rotated, degrees
    :param angle: angle to rotate them by w.r.t origin
    :return:
    """

    xyz = np.copy(pos)
    angle *= (np.pi / 180)  # convert to radians

    R = rotate_z(angle)

    for i in range(np.shape(xyz)[0]):
        xyz[i, :] = np.dot(R, xyz[i, :])

    return xyz


def Rvect2vect(A, B):
    """
    Find rotation of a so that its orientation matches b.
    :param a: vector to be rotated
    :param b: vector to rotate to
    :return: rotation matrix for rotate a to b
    """

    # normalize
    a = A / np.linalg.norm(A)
    b = B / np.linalg.norm(B)

    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = np.dot(a, b)

    v_skew = np.matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    R = np.identity(3) + v_skew + v_skew**2*(1 - c)/(s**2)

    # print(R)
    # cross = np.cross(a, b)  # find vector perpendicular to a and b
    # x = cross / np.linalg.norm(cross)  # normalize
    #
    # y = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    # theta = np.arccos(y)
    #
    # A = np.matrix([[0, -x[2], x[1]], [x[2], 0, -x[0]], [-x[1], x[0], 0]])
    # I = np.identity(3)
    #
    # R = I + np.sin(theta)*A + (1 - np.cos(theta))*(A**2)
    # print(R)
    # exit()

    return np.matrix(R)


def rotate_coords(xyz, R):
    """
    Given a rotation matrix, rotate all points in an array
    :param xyz: xyz coordinates of all positions to be rotated
    :param R: rotate plane
    :return: rotated coordinates
    """
    pos = np.copy(xyz)
    for i in range(np.shape(pos)[0]):
        pos[i, :] = np.dot(R, pos[i, :])

    return pos


def random_orientation(xyz, alignment_vector, placement):
    """
    Randomly orient a water molecule and then place it a desired location
    :param water_xyz: 3d coordinates of molecule
    :param water_alignment_vector: A reference vector to rotate the water molecule about
    :param placement: where to place final water configuration in space
    :return: coordinates of oriented and translated water molecule
    """

    u = np.random.normal(size=3)  # random vector. From normal distribution since sphere
    u /= np.linalg.norm(u)  # normalize

    R = Rvect2vect(alignment_vector, u)  # rotation matrix to align water_alignment_vector with u

    pt = np.random.choice(xyz.shape[0])  # randomly choose reference atom
    xyz -= xyz[pt, :]  # center at origin

    rotated = np.zeros([xyz.shape[0], 3])
    for i in range(xyz.shape[0]):
        rotated[i, :] = np.dot(R, xyz[i, :])

    rotated += placement  # translate to desired location

    return rotated


def rescale(coords, dims):
    """
    rescale coordinates so that cell dimensions are constant over the simulation
    returns rescaled coordinates and average length
    """

    avgdims = np.average(dims, axis=0)
    a = avgdims / dims
    rc = coords

    for it in range(rc.shape[0]):
        for i in range(3):
            rc[it, :, i] *= a[it, i]

    return rc, avgdims


def monoclinic_to_cubic(xyz, theta=60):
    """Convert monoclinic cell to cubic cell"""

    print("transforming coordinates to monoclinic cell (theta={:3.2f} deg)".format(theta*180.0/np.pi))

    coordinates = np.copy(xyz)
    coordinates[..., 1] /= np.sin(theta)
    coordinates[..., 0] -= coordinates[..., 1]*np.cos(theta)

    return coordinates


def fft_3D_monoclinic(xyz, box_vectors, bins, angle=60, sf=False):
    """ Calculate 3D discrete fourier transform for each frame of a trajectory of coordinates in monoclinic unit cell

    :param xyz: frame-by-frames coordinates of atoms whose 3D DFT we want to calculate
    :param box_vectors: matrix of box vectors of shape (nT, 3, 3)
    :param bins: number of bins in each dimension
    :return:
    """

    nT = xyz.shape[0]  # number of frames
    L = np.linalg.norm(box_vectors, axis=2)

    locations, L = rescale(xyz, L)  # make unit cell constant size

    # define histograme bin edges
    x = np.linspace(0, L[0], int(bins[0]))
    y = np.linspace(0, L[1], int(bins[1]))
    z = np.linspace(0, L[2], int(bins[2]))

    zv = [0.0, 0.0, 0.0]  # zero vector

    # put all atoms inside box - works for single frame and multiframe
    for it in range(locations.shape[0]):  # looped to save memory
        locations[it, ...] = np.where(locations[it, ...] < L, locations[it, ...], locations[it, ...] - L)  # get positions in periodic cell
        locations[it, ...] = np.where(locations[it, ...] > zv, locations[it, ...], locations[it, ...] + L)

    # fourier transform loop
    sf = np.zeros([x.size - 1, y.size - 1, z.size - 1])
    for frame in tqdm.tqdm(range(nT), unit=' Frames'):
        H, edges = np.histogramdd(locations[frame, ...], bins=(x, y, z))
        if sf:
            fft = np.fft.fftn(H - H.mean())
            sf += (fft * fft.conjugate()).real
        else:
            sf += np.abs(np.fft.fftn(H)) ** 2

    sf /= nT  # average of all frames

    return sf