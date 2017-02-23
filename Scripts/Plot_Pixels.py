#!/usr/bin/env python

"""
    This script is intended for plotting a 2D image in matplotlib using pixel values from an ASCII format file
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os


def initialize():
    parser = argparse.ArgumentParser(description='Plot a 2D image from an ASCII file')

    parser.add_argument('-f', '--file', default='1periodic_far.asc', help='Name of file to read')
    parser.add_argument('-x', '--xpixel', default=1024, help='Number of pixels in the x direction')
    parser.add_argument('-y', '--ypixel', default=1024, help='Number of pixels in the y direction')
    parser.add_argument('-c', '--cmap', default='Greys', help='Color Scheme for heat map -- more options at '
                                                              'http://matplotlib.org/examples/color/colormaps_reference.html')

    args = parser.parse_args()

    return args


def get_pixels_asc(file, xpixel, ypixel):

    f = open(file, 'r')

    a = []
    for line in f:
        a.append(line)

    f.close()

    pixel_values = np.zeros((int(xpixel), int(ypixel)), dtype=int)

    for i in range(int(ypixel)):
        pixel_values[i, :] = a[i].split()

    return pixel_values


if __name__ == '__main__':
    args = initialize()

    # directory = '/home/bcoscia/Documents/Gromacs/SAXS/saxs_frames_long/pixels'


    # directory = '/home/bcoscia/PycharmProjects/GitHub/Scripts/saxs_frames_far/images'
    #
    # pixel_sum = np.zeros([1024, 1024])
    #
    # count = 0
    # for filename in os.listdir(directory):
    #     if filename.endswith(".asc"):
    #         count += 1
    #         pixel_values = get_pixels('%s/%s' % (directory, filename), '%s' % args.xpixel, '%s' % args.ypixel)
    #         for i in range(1024):
    #             for j in range(1024):
    #                 pixel_values[i, j] -= 40
    #         pixel_sum += pixel_values
    #         continue
    #     else:
    #         continue

    # pixel_values = get_pixels('%s' % args.file, '%s' % args.xpixel, '%s' % args.ypixel)
    # for i in range(1024):
    #     for j in range(1024):
    #         pixel_values[i, j] -= 40
    # plt.figure()
    # im = plt.imshow(pixel_values, cmap='%s' % args.cmap, interpolation='none', vmin=40, vmax=50)

    # import Radial_int_pixels
    #pixels = np.loadtxt(args.file)
    f = open(args.file, 'r')
    a = []
    for line in f:
        a.append(line.split())

    pixels = np.zeros([len(a[0]), len(a)])
    #
    # for i in range(len(a[0])):
    #     print len(a[i])
    #     for j in range(len(a) - 1):
    #         pixels[i, j] = a[i][j]

    for j in range(len(a) - 1):
        for i in range(len(a[0])):
            pixels[i, j] = a[j][i]

    means = []
    mins = []
    for i in range(pixels.shape[1]):
        mins.append(min(pixels[:, i]))
        means.append(np.mean(pixels[:, i]))

    Imin = min(mins)
    print 'Imin = %s' % Imin

    maxes = []
    for i in range(pixels.shape[1]):
        maxes.append(max(pixels[:, i]))

    Imax = max(maxes)
    mean = np.mean(means)
    print 'Imax = %s' % Imax
    print 'mean = %s' % mean

    plt.figure()
    im = plt.imshow(pixels, cmap='%s' % args.cmap, interpolation='none', vmin=1e6, vmax=1e11)
    plt.show()
    # plt.figure()
    # q, theta, intensities = Radial_int_pixels.radial_int(pixel_sum, 1024, .0002, 0.0001, 1.18, 1.54e-10)
    # plt.plot(theta, intensities)
    # plt.show()