#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from past.utils import old_div
import numpy as np
import mdtraj as md
import argparse
from . import tilt
from llclib import file_rw
from scipy import spatial
import matplotlib.pyplot as plt
import math
import os
import tqdm


def initialize():

    parser = argparse.ArgumentParser(description='Measure Membrane Thickness')

    parser.add_argument('-g', '--gro', type=str, default='wiggle.gro', help='Path to input file')
    parser.add_argument('-t', '--traj', type=str, default='traj_whole.xtc', help='Trajectory file (.xtc or .trr)')
    parser.add_argument('-i', '--index', default='tail_index.ndx', type=str, help='Index file containing names of tail carbons')
    parser.add_argument('-c', '--cutoff', default=0.45, type=float, help='Cutoff distance for neighbor search')
    parser.add_argument('--save', action="store_true", help='Save output plots')
    parser.add_argument('--fit', action="store_true", help='Plot a fourier fit')
    parser.add_argument('--start', default=0, type=int, help='start frame')
    parser.add_argument('--end', default=-1, type=int, help='end frame')
    parser.add_argument('--noshow', action="store_true", help='Do not show plots at end')
    parser.add_argument('--suffix', help='Output suffix', default='traj')
    parser.add_argument('--load', action="store_true", help='load previously saved arrays')
    parser.add_argument('--write_gro', action="store_true", help='create .gro file of centroids')

    args = parser.parse_args()

    return args


def tail_centroid(pos, grps):
    """
    :param pos: coordinates of all tail atoms
    :param grps: a list with an entry for each tail. Each entry should have sub-entries for each atom making up the tail
    :return: centroids for each tail
    """

    ngrps = len(grps)
    natoms = pos.shape[1]
    atomsptail = len(grps[0])
    nT = pos.shape[0]
    nmon = natoms / len(grps[0]) / ngrps

    centroids = np.zeros([nT, nmon * ngrps, 3])

    print('Calculating tail centroids')
    for t in tqdm.tqdm(list(range(nT))):
        for i in range(nmon):
            for j in range(ngrps):
                centroid = np.zeros([3])
                for k in range(atomsptail):
                    centroid += pos[t, i*ngrps*atomsptail + j*atomsptail + k, :]
                centroid /= atomsptail
                centroids[t, i*ngrps + j, :] = centroid

    return centroids


def nearest_neighbors(arr, d, lower_limit=0.4):
    """
    :param arr: array of points [npoints, dimension]
    :param d: furthest distance out to search (nm)
    :return: list of list of indices. List i contains all nearest neighbors to position i
    """

    nT = arr.shape[0]
    npts = arr.shape[1]
    nn_list = []
    for t in range(nT):
        nn_list.append([])
        for i in range(npts):
            nn_list[t].append([])  #using lists since this will vary in length

    print('Calculating nearest neighbors')
    for t in tqdm.tqdm(list(range(nT))):
        frame = arr[t, :, :]
        for i in range(npts):
            others = np.delete(frame, i, 0)  # delete self entry or else the nearest neighbor is itself
            nn = spatial.KDTree(others).query(frame[i, :])[1]  # find index of nearest neighbor
            index = np.where(frame == others[nn])  # find the index in others corresponding to nn in arr
            ld = np.linalg.norm(frame[i, :] - frame[index[0][0]])  # calc linear distance between position and nearest neighbor
            while ld <= d:  # make sure its within d and then keep finding nearest neighbors within this range
                if ld >= lower_limit:
                    nn_list[t][i].append(index[0][0])  # add index of nearest neighbor to list i
                others = np.delete(others, nn, 0)  # delete row containing already counted nearest neighbor
                nn = spatial.KDTree(others).query(frame[i, :])[1]  # calculate new nearest neighbor
                index = np.where(frame == others[nn])  # find the index in others corresponding to nn in arr
                ld = np.linalg.norm(frame[i, :] - frame[index[0][0]])  # calc linear distance between position and nearest neighbor

    return nn_list


def angles_ld(arr, indices, normal=[0, 0, 1]):
    """
    :param arr: array of positions
    :param indices: indices of nearest neighbors (output from nearest_neighbors function)
    :return: a distribution of angles and linear distances between nearest neighbors. Angles are w.r.t the plane
    perpendicular to normal. Default is the xy plane
    """

    nT = arr.shape[0]
    npts = arr.shape[1]
    angles = []
    ld = []

    print('Calculating angles between nearest neigbors')
    for t in tqdm.tqdm(list(range(nT))):
        for i in range(npts):
            # for j in range(len(indices[i])):
            for j in indices[t][i]:
                v = arr[t, i, :] - arr[t, j, :]
                vn = np.dot(v, normal)
                nn = np.linalg.norm(normal)
                vv = np.linalg.norm(v)
                angle = np.arcsin(old_div(vn, (nn * vv))) * (old_div(180, np.pi))
                angles.append(angle)
                ld.append(vv)

    return angles, ld


def cn(n, T, bin_locations):
    c = counts*np.exp(-1j*2*n*np.pi*bin_locations/T)
    return old_div(c.sum(),c.size)


def f(x, Nh, T, bin_locations):
    """
    :param x:
    :param Nh:
    :param T: period
    :return:
    """
    f = np.array([2*cn(i, T, bin_locations)*np.exp(1j*2*i*np.pi*x/T) for i in range(1, Nh + 1)])
    return f.sum()


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[args.start:args.end]  # load trajectory

    grps = tilt.read_index(args.index)  # read index file
    ngrps = len(grps)

    all_atoms = []
    for i in range(ngrps):
        for j in range(len(grps[i])):
            all_atoms.append(grps[i][j])

    keep = [a.index for a in t.topology.atoms if a.name in all_atoms]
    tails = t.atom_slice(keep).xyz

    centroids = tail_centroid(tails, grps)

    if args.write_gro:
        file_rw.write_gro_pos(centroids[-1, :, :], 'centroids.gro', name='NA')

    if not args.load:
        nlist = nearest_neighbors(centroids, args.cutoff)
    else:
        np.load('nlist.npy')

    if args.save and not args.load:
        np.save('nlist', nlist)

    angles, ld = angles_ld(centroids, nlist)

    np.savez_compressed('angles_ld_%s.npz' % args.suffix, angles=angles, ld=ld, nlist=nlist)

    angles = [value for value in angles if not math.isnan(value)]
    plt.figure(1)
    nbins = 45
    (counts, bins, patches) = plt.hist(angles, bins=nbins)

    if args.fit:
        bin_locations = np.linspace(-90, 90, nbins)
        fapprox = np.array([f(t, 8, 180, bin_locations).real for t in counts]) + np.mean(bins)
        plt.plot(bin_locations, fapprox)

    if args.save:
        plt.savefig('angles_%s.png' % args.suffix)

    plt.figure(2)
    plt.hist(ld)
    if args.save:
        plt.savefig('distances_%s.png' % args.suffix)

    if not args.noshow:
        plt.show()