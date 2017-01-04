#!/usr/bin/python

import argparse
import Radial_int_pixels
import numpy as np
import matplotlib.pyplot as plt
import Atom_props
import time


def initialize():
    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-f', '--file', default='10periodic.asc', help='Name of file to read')
    parser.add_argument('-p', '--pixels', default=1024, help='Number of pixels in each direction (i.e. p x p)')
    parser.add_argument('-r', '--dr', default=0.0002, help='Step size for integration')
    parser.add_argument('-d', '--dim', default=0.0001, help='Dimension of pixel (assuming a square pixel) [m]')
    parser.add_argument('-w', '--wavelength', default=1.54, help='Wavelength of X-rays, [m]')
    parser.add_argument('-D', '--dist', default=0.1, help='Distance from sample center to detector center')
    parser.add_argument('-i', '--incident', default=[0, 0, 1], help='Incident x_ray vector - input as list')

    args = parser.parse_args()
    return args


def detector_coords(pixels, dim, dist):
    """

    :param pixels: The number of pixels in the x and y dimension - assumed to be a square detector
    :param dim: Pixel dimensions (m)
    :return: A: The x, y, z positions of each pixel on the detector relative to the origin placed at the sample location
    """
    pixels = int(pixels)
    dim = float(dim)
    dist = float(dist)

    A = np.zeros([3, pixels, pixels])

    m = np.shape(A)[1]
    n = np.shape(A)[2]
    shift = m - 1

    for i in range(m):
        for j in range(n):
            A[0, i, j] = (2 * j - shift) * dim / 2
            A[1, i, j] = (- 2 * i + shift) * dim / 2

    A[2, :, :] = dist
    return A


def k_vectors(detector, wavelength):
    """
    :param: The x, y, z positions of each pixel on the detector relative to the origin placed at the sample location
    :return: k_vect: The normalized vectors pointing from the sample center to each pixel
    """

    origin = np.zeros([3])
    pixels = np.shape(detector)[1]
    K = 2 * np.pi / float(wavelength)  # desired magnitude of the wave-vector

    k_vect = np.zeros([3, pixels, pixels])
    for i in range(pixels):
        for j in range(pixels):
            k_vect[:, i, j] = K * detector[:, i, j] / np.linalg.norm(detector[:, i, j])

    return k_vect

if __name__ == '__main__':
    args = initialize()
    pixels = int(args.pixels)

    detector_pts = detector_coords(pixels, '%s' % args.dim, '%s' % args.dist)

    k_vect = k_vectors(detector_pts, '%s' % args.wavelength)

    # pos = np.load('NApos_array612ns')
    # id = np.load('identity_array612ns')

    # pos = pos[:, :, 0]

    pos = np.load('NA_positions_5_images')
    atoms = np.shape(pos)[1]
    # print atoms
    incident = np.array(args.incident)

    # Switch y and z coordinate
#    for i in range(atoms):
 #       temp = pos[1, i]
  #      pos[1, i] = pos[2, i]
  #      pos[2, i] = temp

    Intensities = np.zeros([pixels, pixels])

    for i in range(pixels):
        print '{:2.1f} percent done'.format(100.0 * (i + 1) / pixels)
        if i != 0:
            tm = (end - start) * (pixels - i - 1)
            hours = int(np.floor(tm / 3600))
            minute = int(np.floor((tm - 3600 * hours) / 60))
            seconds = int(tm - (3600 * hours) - (60 * minute))
            print 'Approximate Time remaining: %02d:%02d:%02d' % (hours, minute, seconds)
        start = time.time()
        for j in range(pixels):
            Re = 0
            Im = 0
            for k in range(atoms):
                # Z = Atom_props.mass[id[0, k]]
                Z = 22.98977
                rq = np.dot(pos[:, k], (k_vect[:, i, j] - incident))
                Re += Z*np.cos(rq)
                Im += Z*np.sin(rq)
            Intensities[i, j] += Re ** 2 + Im ** 2
        end = time.time()

    mins = []
    for i in range(1024):
        mins.append(min(Intensities[:, i]))

    Imin = min(mins)

    maxes = []
    for i in range(pixels):
        maxes.append(max(Intensities[:, i]))

    Imax = max(maxes)
    print 'Imax = %s' % Imax

    #im = plt.imshow(Intensities, cmap='Greys', interpolation='none', vmin=Imin, vmax=Imax)

    f = open('Intensities', 'w')
    np.save(f, Intensities)
    f.close()

    plt.show()
