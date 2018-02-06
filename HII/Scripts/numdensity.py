#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import math
from pymbar import timeseries


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the number density of selected groups along axis')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='wiggle.trr', help='Trajectory file (.trr or .xtc)')
    parser.add_argument('-d', '--axis', default='z', help='Axis along which to calculate number density. If you put '
                                                          'anything other than x, y or z, it will default to z')
    parser.add_argument('-a', '--atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], type=str, help='List of atoms of interest')
    parser.add_argument('-c', '--center', default='yes', help='Set this to yes if you want to calculate the number'
                                                              'density based on the centers of the selected atoms')
    parser.add_argument('-b', '--bin', default=.1, type=float, help='bin size (nm)')
    parser.add_argument('--fourier', action="store_true", help='Create a power spectrum')
    parser.add_argument('--begin', default=0, type=int, help='start frame')
    parser.add_argument('--end', default=-1, type=int, help='end frame')
    parser.add_argument('--entropy', action="store_true")
    parser.add_argument('--noshow', action="store_true")
    parser.add_argument('--trajectory', action="store_false")
    parser.add_argument('--kb', action="store_true", help='Calculate kullback-leibler divergence')
    parser.add_argument('--layers', action="store_true", help='Calculate the number of layers by binning')
    parser.add_argument('-bb', '--binbins', type=int, help='number of bins')
    parser.add_argument('--shift', default=0, type=int, help='shift starting point for binning to the right')

    args = parser.parse_args()

    return args


def centers(pos, atoms):
    """
    Find the average coordinates based on coordinates of selected atoms. (useful for rings i.e. benzene)
    :param pos: numpy array with xyz coordinates of selected atoms for all frames [frames, no_atoms, xyz coords]
    :param atoms: list of selected atoms
    :return: the coordinates for the centers of all selected atoms
    """

    if len(pos.shape) == 3:
        frames = pos.shape[0]
        natoms = pos.shape[1]
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
    else:  # single frame case
        natoms = pos.shape[0]
        ncenters = natoms / len(atoms)
        frames = 1

        c = np.zeros([1, ncenters, 3])
        nselect = len(atoms)

        for i in range(frames):
            for j in range(ncenters):
                sum = np.zeros([3])
                for k in range(nselect):
                    sum += pos[j*nselect + k, :]
                sum /= nselect
                c[i, j, :] = sum

    return c


def d_dist(pos, box, bin):
    """
    Find the distribution of distances of components from each (basically a radial distribution function)
    :param pos: A numpy array of positions of components of interest
    :param box: A numpy array of box vectors used to find the maximum distance
    :return: Counts of components at a distance vs. distance
    """
    nT = box.shape[0]
    nV = box.shape[1]
    nA = pos.shape[1]
    max_d = 0
    for i in range(nT):
        for j in range(nV):
           if box[i, j] > max_d:
               max_d = box[i, j]

    r = np.linspace(0, max_d, max_d/bin + 1)
    dist = np.zeros(len(r))
    bins = len(dist)

    for i in range(nT):
        for j in range(nA):
            for k in range(nA):
                if j != k:
                    d = np.linalg.norm(pos[i, j, :] - pos[i, k, :])
                    bin = int(bins * (d / max_d))  # rounds down
                    dist[bin] += 1

    return dist, r


def density(pos, axis, bin, box, sum=True, smooth_factor=1):
    """
    Calculate the average number density of components along an axis
    :param pos: a numpy array with xyz coordinates of selected atoms for all frames
    :param axis: which axis to split up (x, y or z)
    :param bin: size of the bins which the axis will be split into
    :param box: the box
    :param sum: if sum = 'yes' then all frames will be added together and averaged
    :param smooth_factor = x: record measurements every x frames (only relevant for trajectories)
    :return: density as a function of distance along axis
    """

    nT = pos.shape[0]  # number of trajectory points
    nA = pos.shape[1]  # number of atoms

    if nT == 1:
        b = box[axis] / 2
    else:
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
    L = len(x)
    if sum:
        d = np.zeros([L])
    else:
        d = np.zeros([int(nT/smooth_factor), L])
        xT = np.zeros([int(nT/smooth_factor), L])

    for i in range(int(nT/smooth_factor)):

        if nT != 1:
            b = box[i*smooth_factor, axis] / 2  # center of the membrane with respect to axis
        else:
            b = box[axis] / 2

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

        if not sum:
            xT[i, :min(L, len(x))] = x[:min(L, len(x))]

        for j in range(nA):
            a = pos[i*smooth_factor, j, axis]

            bin_no = int(len(x) / 2 + np.floor((a - b)/bin))

            if sum:
                if bin_no < L:
                    d[bin_no] += 1
            else:
                if bin_no < L:
                    d[i, bin_no] += 1

    if sum:
        for i in range(len(d)):  # take the average
            d[i] /= (nT / smooth_factor)

    if sum:
        return x, d
    else:
        return xT, d


def power_spectrum(data, bin):
    """
    Compute the power spectrum of the data (find dominant frequencies in the fourier series fitting discrete data)
    :param data: Data to be fourier transformed (1D numpy array)
    :param bin: bin size (nm)
    :return: power spectrum of data (ps) with corresponding frequencies (freqs), the max frequency (max) and
    """

    x = np.linspace(0, bin*len(data), len(data))

    # subtract the mean from all data points
    avg = np.mean(data)
    data = np.array([abs(i - avg) for i in data])

    N = data.size

    data = data[int(.1*N):int(.9*N)]
    data = data - np.mean(data)
    ps = np.abs(np.fft.fft(data))**2

    freqs = np.fft.fftfreq(data.size)
    idx = np.argsort(freqs)

    max_freq = np.argmax(np.abs(np.fft.fft(data)))
    freq = freqs[max_freq]

    # modify things so they'll plot nicely and in the correct units
    freqs = freqs[idx] / bin
    ps = ps[idx]
    max = abs(freq / bin)  # maximum frequency in hertz

    return ps, freqs, max


def entropy(d, traj=False):
    """
    Calculate the entropy of a discrete data series using information theory
    :param x: possible values
    :param d: distribution
    """

    if traj:  # return time depenedent entropy based on a trajectory

        nT = d.shape[0]  # number of points in trajectory
        nbins = d.shape[1]
        H = np.zeros([nT])

        for t in range(nT):
            # normalize
            frame = np.zeros([nbins])
            frame = d[t, :]
            frame /= sum(frame)
            frame = np.trim_zeros(frame, 'fb')

            for b in range(len(frame)):
                if frame[b] != 0:  # https://stats.stackexchange.com/questions/57069/alternative-to-shannons-entropy-when-probability-equal-to-zero
                    H[t] += -frame[b]*math.log(frame[b])

    else:
        # normalize the distribution
        d /= sum(d)
        d = np.trim_zeros(d, 'fb')  # get rid of zeros at front and back of distribution

        H = 0  # entropy
        for i in range(d.shape[0]):
            if d[i] != 0:
                H += -d[i]*(math.log(d[i]))

    return H


def kullback_leiblier(p, q):
    """
    Find the kullback-leiblier divergence from q to p

    :param p: probability distribution 1
    :param q: probability distribution 2
    :return: kb : kullback-leiblier divergence - the amount of information lost when q is used to approximate p
    """
    # normalize
    p /= sum(p)
    q /= sum(q)

    nbins = p.shape[0]

    kb = 0
    for i in range(nbins):
        if p[i] != 0 and q[i] != 0:
            kb += p[i]*math.log(p[i]/q[i])

    return kb


if __name__ == '__main__':

    args = initialize()

    t = md.load('%s' % args.traj, top='%s' % args.gro)[args.begin:args.end]

    atoms = args.atoms
    keep = [a.index for a in t.topology.atoms if a.name in atoms]  # restrict trajectory to chosen atoms
    pos = t.atom_slice(keep).xyz  # get just the coordinates
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

    # uncomment to see radial distribution
    # dist, r = d_dist(c, box, args.bin)  # this takes forever. Use a gromacs tool. Probably gmx rdf
    # plt.bar(r, dist, args.bin)
    # plt.show()

    x, d = density(c, axis, args.bin, box, sum=args.trajectory)
    nfig = 1

    if args.layers:
        # create a number of bins equal to the number of layers we are guessing
        # leading_zeros = 0
        # while d[leading_zeros] == 0:
        #     leading_zeros += 1
        #
        # trailing_zeros = -1
        # while d[trailing_zeros] == 0:
        #     trailing_zeros -= 1
        #
        # trailing_zeros += d.shape[0]
        # nnonzero = trailing_zeros - leading_zeros

        shift = args.shift
        density = np.trim_zeros(d, 'fb')  # get ride of leading and trailing zeros of density distribution
        nvalues = density.shape[0]
        values_per_bin = nvalues / args.binbins  # number of values in each bin
        bins = np.zeros([values_per_bin])
        for i in range(args.binbins):
            for j in range(values_per_bin):
                bins[j] += density[(i*values_per_bin + j + shift) % nvalues]

        width = args.bin * values_per_bin
        bin_width = x[1] - x[0]
        x_values = np.linspace(0, width, values_per_bin)
        plt.bar(x_values, bins, bin_width)
        plt.show()
        exit()

    # to avoid confusion reading this later
    traj = True
    if len(x.shape) == 1:
        traj = False

    if traj:
        bins = d.shape[1]
    else:
        bins = d.shape[0]  # because trajectory

    if args.fourier:

        ps, freqs, max_freq = power_spectrum(d, args.bin)

        print('Maximum frequency: %s cycles/nm' % max_freq)

        plt.figure(nfig)
        plt.plot(freqs, ps)
        plt.suptitle('Power Spectrum', fontsize=16)
        plt.title('Bin size = %s nm' % args.bin, fontsize=12)
        plt.xlabel('Frequency')
        plt.ylabel('Fourier transformed data squared')

        nfig += 1

    if args.entropy:

        H_real = entropy(d, traj)

        H_uniform = entropy(np.ones(bins), traj=False)

        if traj:
            H_diff = np.zeros(H_real.shape)
            for frame in range(H_real.shape[0]):
                H_diff[frame] = H_uniform - H_real[frame]

            H_diff /= max(H_diff)
            equil = timeseries.detectEquilibration(H_diff)[0]
            print("Equilibration after %d ns" % (t.time[equil] / 1000))
            avg_order = np.mean(H_diff[equil:])

        else:
            print('Entropy difference from uniform distribution : %s' % (H_uniform - H_real))

    if args.kb:
        # Calculate Kullback-Leibler divergence
        if traj:
            x_kb, d_kb = density(c, axis, args.bin, box, sum=True)
            kb_fwd = kullback_leiblier(np.ones(bins), d_kb)
            kb_bwd = kullback_leiblier(d_kb, np.ones(bins))

        else:
            kb_fwd = kullback_leiblier(np.ones(bins) / bins, d)
            kb_bwd = kullback_leiblier(d, np.ones(bins) / bins)

        print('kb_fwd = %.3f' % kb_fwd)
        print('kb_bwd = %.3f' % kb_bwd)

    if not args.noshow:

        if traj and args.entropy:
            plt.figure(nfig)
            plt.plot(t.time, H_diff)
            # plt.title('Entropy vs. time')
            # plt.ylabel('$/delta$ S between real and uniform distribution')
            plt.title('Equilibrium order fraction : %.2f' % avg_order)
            plt.ylabel('Order Fraction')
            plt.xlabel('Time (ps)')
            plt.ylim(ymin=0)
            nfig += 1

        bin_width = max(x[:d.shape[0]]) / d.shape[0]
        plt.figure(nfig)
        plt.bar(x[:d.shape[0]], d, bin_width)
        plt.suptitle('Line number density of components along %s axis' % args.axis, fontsize=16)
        plt.title('Bin size = %s nm' % args.bin, fontsize=12)
        plt.xlabel('Distance into membrane, z direction (nm)')
        plt.ylabel('NA count')
        plt.show()