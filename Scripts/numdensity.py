#!/usr/bin/python

import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import math


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the number density of selected groups along axis')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='wiggle.trr', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-d', '--axis', default='z', help='Axis along which to calculate number density. If you put '
                                                          'anything other than x, y or z, it will default to z')
    parser.add_argument('-a', '--atoms', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='List of atoms of interest')
    parser.add_argument('-c', '--center', default='yes', help='Set this to yes if you want to calculate the number'
                                                              'density based on the centers of the selected atoms')
    parser.add_argument('-b', '--bin', default=.4, type=float, help='bin size (nm)')

    args = parser.parse_args()

    return args


def centers(pos, atoms):
    """
    Find the average coordinates based on coordinates of selected atoms. (useful for rings i.e. benzene)
    :param pos: numpy array with xyz coordinates of selected atoms for all frames
    :param atoms: list of selected atoms
    :return: the coordinates for the centers of all selected atoms
    """
    frames = np.shape(pos)[0]
    natoms = np.shape(pos)[1]
    ncenters = natoms / len(atoms)  # the number of centers that will be calculated
    c = np.zeros([frames, ncenters, 3])
    nselect = len(atoms)  # the number of selected atoms

    for i in range(frames):
        for j in range(ncenters):
            sum = np.zeros([3])
            for k in range(nselect):
                sum += pos[i, j*nselect + k, :]
            sum /= nselect
            c[i, j, :] = sum

    return c


def density(pos, axis, bin, box):
    """
    Calculate the number density of components along an axis (last frame only for now)
    :param pos: a numpy array with xyz coordinates of selected atoms for all frames
    :param axis: which axis to split up (x, y or z)
    :param bin: size of the bins which the axis will be split into
    :param box: the box
    :return: density as a function of distance along axis
    """

    b = box[-1, axis]
    print b
    x = np.linspace(0, b, b/bin + 1)
    d = np.zeros([len(x)])
    for i in range(np.shape(pos)[1]):
        bin_no = np.floor(pos[-1, i, axis] / bin)
        d[bin_no] += 1

    return x, d

if __name__ == '__main__':

    args = initialize()

    # load trajectory
    t = md.load('%s' % args.traj, top='%s' % args.gro)
    atoms = args.atoms
    keep = [a.index for a in t.topology.atoms if a.name in atoms]  # restrict trajectory to chosen atoms
    t.restrict_atoms(keep)
    pos = t.xyz  # get just the coordinates
    box = t.unitcell_lengths  # get the unit cell lengths

    frames = np.shape(pos)[0]

    # find out along which axis we are going to analyze
    if args.axis == 'x':
        axis = 0
    elif args.axis == 'y':
        axis = 1
    else:
        axis = 2

    c = centers(pos, args.atoms)
    x, d = density(c, axis, float(args.bin), box)


    d = np.trim_zeros(d)
    x = np.linspace(0, args.bin*len(d), len(d))
    f = open('d', 'w')
    np.save(f, d)
    f.close()
    f = open('x', 'w')
    np.save(f, x)
    f.close()
    # avg = np.mean(d)
    # d = np.array([abs(i - avg) for i in d])
    fft = np.fft.fft(x)

    # N = fft.shape[-1]
    # A = np.zeros([N])
    # phi = np.zeros([N])
    # f = np.zeros([N])
    # print np.fft.fftfreq(fft.shape, d=0.1)
    # for i in range(N):
    #     A[i] = np.sqrt(fft[i].real ** 2 + fft[i].imag ** 2)
    #     phi[i] = math.atan2(fft[i].imag, fft[i].real)
    #     f[i] = A[i]/(2*N)

    # fx = np.zeros([N])
    # for i in range(N):
    #     for j in range(N):
    #         fx[i] += A[]

    N = d.size

    data = d[int(.1*N):int(.9*N)]
    data = data - np.mean(data)
    x = x[int(.1*N):int(.9*N)]
    ps = np.abs(np.fft.fft(data))**2

    time_step = args.bin
    freqs = np.fft.fftfreq(data.size, time_step)
    idx = np.argsort(freqs)

    plt.plot(freqs[idx], ps[idx])
    # plt.plot(x, np.fft.fft(data))
    plt.title('Power Spectrum')
    plt.xlabel('Frequency')
    plt.ylabel('Fourier transformed data squared')
    plt.show()

    # plt.plot(x, d)
    # plt.xlabel('Distance into membrane, z direction (nm)')
    # plt.ylabel('Count of benzene centers')
    # plt.show()