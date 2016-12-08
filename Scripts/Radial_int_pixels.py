#!/usr/bin/python

"""
    This script is intended to radially integrate a 2D image defined by a numpy array of pixel intensities
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import Plot_Pixels


def initialize():

    parser = argparse.ArgumentParser(description='Radially integrate a 2D image defined by a numpy array of pixel '
                                                 'intensities')

    parser.add_argument('-f', '--file', default='10periodic.asc', help='Name of file to read')
    parser.add_argument('-p', '--pixels', default=1024, help='Number of pixels in each direction (i.e. p x p)')
    parser.add_argument('-r', '--dr', default=0.0002, help='Step size for integration')
    parser.add_argument('-d', '--dim', default=0.0001, help='Dimension of pixel (assuming a square pixel) [m]')
    parser.add_argument('-w', '--wavelength', default=1e-10, help='Wavelength of X-rays, [m]')
    parser.add_argument('-D', '--dist', default=0.1, help='Distance from sample center to detector center')

    args = parser.parse_args()

    return args


def dist_mat(pixels, dim):

    dim = float(dim)
    pixels = int(pixels)

    A = np.zeros([2, pixels, pixels])

    m = np.shape(A)[1]
    n = np.shape(A)[2]
    shift = m - 1

    for i in range(m):
        for j in range(n):
            A[0, i, j] = (2 * j - shift) * dim / 2
            A[1, i, j] = (- 2 * i + shift) * dim / 2

    distance_matrix = np.zeros([pixels, pixels])
    for i in range(m):
        for j in range(n):
            x = A[:, i, j]
            dist = np.sqrt(np.dot(x, x))
            distance_matrix[i, j] = dist

    maxes = []
    for i in range(m):
        maxes.append(max(distance_matrix[:, i]))

    max_dist = max(maxes)

    return distance_matrix, max_dist


def angle_mat(distance_matrix, dist, wavelength):

    m = np.shape(distance_matrix)[0]
    n = np.shape(distance_matrix)[1]
    d = float(dist)

    angles = np.zeros([m, n])

    for i in range(m):
        for j in range(n):
            p = distance_matrix[i, j]
            angles[i, j] = np.arctan(p / d)

    maxes = []
    for i in range(m):
        maxes.append(max(angles[:, i]))

    theta_max = max(maxes)
    print theta_max
    d_mat = np.zeros([m, n])
    w = float(wavelength)
    m2A_conv = 1e10

    for i in range(m):
        for j in range(n):
            d_mat[i, j] = w / (2 * np.sin(angles[i, j])) * m2A_conv

    maxes = []
    for i in range(m):
        maxes.append(max(d_mat[:, i]))

    dmax = np.ceil(max(maxes))

    q_mat = np.zeros([m, n])

    for i in range(m):
        for j in range(n):
            # q_mat[i, j] = 2*np.pi / d_mat[i, j]
            q_mat[i, j] = 4*np.pi*np.sin(angles[i, j] / 2) / (w * m2A_conv)

    maxes = []
    for i in range(m):
        maxes.append(max(q_mat[:, i]))

    qmax = np.ceil(max(maxes))

    return angles, theta_max, q_mat, qmax, d_mat, dmax


def radial_int(pixel_intensities, pixels, dr, dim, dist, wavelength):

    dr = float(dr)
    pixels = int(pixels)
    dim = float(dim)
    dist = float(dist)

    distance_mat, max_dist = dist_mat(pixels, dim)
    print max_dist
    angles, theta_max, q_mat, qmax, d_mat, dmax = angle_mat(distance_mat, dist, wavelength)

    intensities = np.zeros([int(max_dist/dr)])
    bins = len(intensities)
    print bins

    x = np.linspace(0, max_dist, bins)
    qs = np.linspace(0, qmax, bins)
    thetas = np.linspace(0, theta_max, bins)

    for i in range(pixels):
        for j in range(pixels):
            intensity = pixel_intensities[i, j]
            distance = distance_mat[j, i]  # indices flipped in this matrix ...
            bin_no = np.floor(distance/max_dist * bins)
            if bin_no == bins:
                intensities[bin_no - 1] += intensity
            else:
                intensities[bin_no] += intensity

    # for i in range(pixels):
    #     for j in range(pixels):
    #         intensity = pixel_intensities[i, j]
    #         q = q_mat[i, j]  # indices flipped in this matrix ...
    #         bin_no = np.floor(q/qmax * bins)
    #         intensities[bin_no] += intensity

    return qs, thetas, intensities


if __name__ == '__main__':
    args = initialize()
    pixel_intensities = Plot_Pixels.get_pixels('%s' % args.file, '%s' % args.pixels, '%s' % args.pixels)
    q, theta, intensity = radial_int(pixel_intensities, '%s' % args.pixels, '%s' % args.dr, '%s' % args.dim,
                                     '%s' % args.dist, '%s' % args.wavelength)
    plt.plot(q, intensity)
    plt.show()