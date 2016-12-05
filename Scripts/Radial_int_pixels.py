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
    parser.add_argument('-r', '--dr', default=.5, help='Step size for integration')

    args = parser.parse_args()

    return args


def dist_mat(pixels):

    pixels = int(pixels)

    A = np.zeros([2, pixels, pixels])

    m = np.shape(A)[1]
    n = np.shape(A)[2]
    shift = m - 1

    for i in range(m):
        for j in range(n):
            A[0, i, j] = 2 * j - shift
            A[1, i, j] = - 2 * i + shift

    distance_matrix = np.zeros([pixels, pixels])
    for i in range(m):
        for j in range(n):
            x = A[:, i, j]
            dist = np.sqrt(np.dot(x, x))
            distance_matrix[i, j] = dist

    maxes = []
    for i in range(m):
        maxes.append(max(distance_matrix[:, i]))

    max_dist = np.ceil(max(maxes))

    return distance_matrix, max_dist


def radial_int(pixel_intensities, pixels, dr):

    dr = float(dr)
    pixels = int(pixels)

    d_mat, max_dist = dist_mat(pixels)

    intensities = np.zeros([int(max_dist/dr)])
    bins = len(intensities)
    x = np.linspace(0, max_dist, bins)

    for i in range(pixels):
        for j in range(pixels):
            intensity = pixel_intensities[i, j]
            distance = d_mat[j, i]  # indices flipped in this matrix ...
            bin_no = np.floor(distance/max_dist * bins)
            intensities[bin_no] += intensity

    return x, intensities


if __name__ == '__main__':
    args = initialize()
    pixel_intensities = Plot_Pixels.get_pixels('%s' % args.file, '%s' % args.pixels, '%s' % args.pixels)
    x, intensities = radial_int(pixel_intensities, '%s' % args.pixels, '%s' % args.dr)
    plt.plot(x, intensities)
    plt.show()