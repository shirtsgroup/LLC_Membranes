#!/usr/bin/python

"""
    This script is intended for plotting a 2D image in matplotlib using pixel values from an ASCII format file
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import os


def initialize():
    parser = argparse.ArgumentParser(description='Plot a 2D image from an ASCII file')

    parser.add_argument('-f', '--file', default='10periodic.asc', help='Name of file to read')
    parser.add_argument('-x', '--xpixel', default=1024, help='Number of pixels in the x direction')
    parser.add_argument('-y', '--ypixel', default=1024, help='Number of pixels in the y direction')
    parser.add_argument('-c', '--cmap', default='Greys', help='Color Scheme for heat map -- more options at '
                                                              'http://matplotlib.org/examples/color/colormaps_reference.html')

    args = parser.parse_args()

    return args


def get_pixels(file, xpixel, ypixel):

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

    directory = '/home/bcoscia/Documents/Gromacs/SAXS/saxs_frames_long/pixels'

    pixel_sum = np.zeros([1024, 1024])

    for filename in os.listdir(directory):
        if filename.endswith(".asc"):
            pixel_values = get_pixels('%s' % args.file, '%s' % args.xpixel, '%s' % args.ypixel)
            for i in range(1024):
                for j in range(1024):
                    pixel_values[i, j] -= 40
            pixel_sum += pixel_values
            continue
        else:
            continue

    # pixel_values = get_pixels('%s' % args.file, '%s' % args.xpixel, '%s' % args.ypixel)
    # for i in range(1024):
    #     for j in range(1024):
    #         pixel_values[i, j] -= 40
    plt.figure()
    im = plt.imshow(pixel_sum, cmap='%s' % args.cmap, interpolation='none', vmin=0, vmax=10)

    import Radial_int_pixels

    plt.figure()
    x, intensities = Radial_int_pixels.radial_int(pixel_sum, 1024, 2)
    print x
    plt.plot(x, intensities)
    plt.show()