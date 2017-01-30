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
    parser.add_argument('-a', '--atoms', default=['C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22',
                 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37',
                 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48'], help='List of atoms of interest')
    parser.add_argument('-c', '--center', default='yes', help='Set this to yes if you want to calculate the number'
                                                              'density based on the centers of the selected atoms')
    parser.add_argument('-b', '--bin', default=.01, type=float, help='bin size (nm)')

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
    Calculate the average number density of components along an axis
    :param pos: a numpy array with xyz coordinates of selected atoms for all frames
    :param axis: which axis to split up (x, y or z)
    :param bin: size of the bins which the axis will be split into
    :param box: the box
    :return: density as a function of distance along axis
    """

    nT = pos.shape[0]  # number of trajectory points
    nA = pos.shape[1]  # number of atoms

    b = box[0, axis] / 2  # middle of the box with respect to axis
    count = 0
    x = []
    while (b - bin*count) >= 0:  # increment from box center to bottom of box using bin size. (POTENTIAL ERROR here)
        x.append(b - bin*count)
        count += 1
    # potential future errors on the bounds of the box. It is assumed that the bottom of the box is at 0 and the top is
    # at whatever the box[i, axis] is equal to
    count = 1
    while (b + bin*count) <= 2*b:  # increment from box center to bottom of box using bin size. (POTENTIAL ERROR here)
        x.append(b + bin*count)
        count += 1

    x.sort()  # sort the list in order
    x = np.array(x)
    d = np.zeros([len(x) - 1])

    for i in range(nT):

        b = box[i, axis] / 2  # center of the membrane with respect to axis

        count = 0
        x = []
        while (b - bin*count) >= 0:
            x.append(b - bin*count)
            count += 1

        count = 1
        while (b + bin*count) <= 2*b:
            x.append(b + bin*count)
            count += 1

        x.sort()
        x = np.array(x)

        for j in range(nA):
            a = pos[i, j, axis]
            print a
            bin_no = len(x) / 2 + np.floor((a - b)/bin)
            d[bin_no] += 1

    for i in range(len(d)):  # take the average
        d[i] /= nT

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
    print np.shape(box)
    exit()
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

    N = d.size

    data = d[int(.1*N):int(.9*N)]
    data = data - np.mean(data)
    ps = np.abs(np.fft.fft(data))**2

    time_step = args.bin
    freqs = np.fft.fftfreq(data.size)
    idx = np.argsort(freqs)

    max_freq = np.argmax(np.abs(np.fft.fft(data)))
    freq = freqs[max_freq]

    freq_in_hertz = abs(freq / args.bin)

    print 'Maximum frequency: %s cycles/nm' % freq_in_hertz

    plt.figure(1)
    plt.plot(freqs[idx] / args.bin, ps[idx])
    plt.suptitle('Power Spectrum', fontsize=16)
    plt.title('Bin size = %s nm' % args.bin, fontsize=12)
    plt.xlabel('Frequency')
    plt.ylabel('Fourier transformed data squared')

    plt.figure(2)
    plt.plot(x, d)
    plt.suptitle('Line number density of components along %s axis' % args.axis, fontsize=16)
    plt.title('Bin size = %s nm' % args.bin, fontsize=12)
    plt.xlabel('Distance into membrane, z direction (nm)')
    plt.ylabel('NA count')
    plt.show()