#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import time
import matplotlib.pyplot as plt
import pymbar
import tqdm


def initialize():

    parser = argparse.ArgumentParser(description='Calculate nematic order parameter')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file. Make sure to '
                                                                            'preprocess with gmx trjconv -pbc whole')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate file')
    parser.add_argument('-p', '--plane', default=['C', 'C2', 'C4'], help='Names of atoms making up planes whose normal'
                                                                         'we are interested in')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Start frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='End frame')
    parser.add_argument('-bins', default=50, type=int, help='Number of bins in final histogram')
    parser.add_argument('-blocks', type=int, help='Number of blocks for block averaging')
    parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials')

    args = parser.parse_args()

    return args


def normals(t, plane):
    """
    Calculate normal vectors to planes defined by atoms in 'plane'
    :param t: mdtraj trajectory object
    :param plane: names of atoms in plane
    :return: normal vectors to all planes for each frame
    """

    keep = [a.index for a in t.topology.atoms if a.name in args.plane]
    pos = t.xyz[:, keep, :]
    nplanes = pos.shape[1] // len(args.plane)
    nT = pos.shape[0]

    n = np.zeros([nT, nplanes, 3])

    # start = time.time()
    for t in range(nT):
        v1 = pos[t, 1::3, :] - pos[t, ::3, :]
        v2 = pos[t, 2::3, :] - pos[t, 1::3, :]
        n[t, :, :] = np.cross(v1, v2)
    # print(time.time() - start)

    # for a future demonstration. This is 160x slower than the above loop. Can show this whole script up to here as an
    # example. I.e. use atom_slice to show how slow that is
    # n_slow = np.zeros([nT, nplanes, 3])
    # start = time.time()
    # for t in range(nT):
    #     for p in range(nplanes):
    #         v1 = pos[t, p*3 + 1, :] - pos[t, p*3, :]
    #         v2 = pos[t, p*3 + 2, :] - pos[t, p*3 + 1, :]
    #         n_slow[t, p, :] = np.cross(v1, v2)
    # print(time.time() - start)

    return n


def angles(vector, plane=None):
    """
    Calculate angle between a vector and plane
    :param vector: vector(s) whose angle with 'plane' we want
    :param plane: vector normal to plane
    :return:
    """

    nT = vector.shape[0]
    costheta = np.zeros_like(vector[:, :, 0])
    pn = np.linalg.norm(plane)  # norm of normal vector to plane
    lim = np.cos(np.pi / 2)
    # print(lim)
    # exit()

    # start = time.time()
    for i in range(nT):
        costheta[i, :] = np.dot(vector[i, :, :], plane) / (np.linalg.norm(vector[i, :, :], axis=1)*pn)
        # for j in tqdm.tqdm(range(costheta.shape[1])):  # for symmetric molecules (can turn over and molecule is the same)
        #     if lim < costheta[i, j]:
        #         costheta[i, j] -= 1
        #     elif -lim > costheta[i, j]:
        #         costheta[i, j] += 1

    # print(time.time() - start)

    # slower loop
    # costheta2 = np.zeros_like(vector[:, :, 0])
    # nplanes = vector.shape[1]
    # start = time.time()
    # for i in range(nT):
    #     for j in range(nplanes):
    #         costheta2[i, j] = np.dot(vector[i, j, :], plane) / (np.linalg.norm(vector[i, j, :])*pn)
    # print(time.time() - start)

    return costheta


def nematic(costheta):
    """
    Calculate the nematic order parameter of a liquid crystal system
    See p. 168 of : http://www.dsf.unica.it/~fiore/libricorsoptr/Chaikin%20Lubensky%20-%20Principles%20of%20Condensed%20Matter%20Physics.jb2.pdf
    :param costheta: cosine of the angle made with director vector (get from "angles" function above)
    :return: nematic order parameter at each frame
    """

    S = np.zeros([costheta.shape[0]])  # calculate S for each frame
    for i in range(costheta.shape[0]):
        S[i] = 0.5*np.sum(3*a[i, :]**2 - 1)

    S /= costheta.shape[1]

    return S


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[args.begin:args.end]
    nT = t.n_frames

    n = normals(t, args.plane)

    a = angles(n, [0, 0, 1])

    plt.hist(a.flatten(), bins=100)
    plt.show()
    exit()

    S = nematic(a)
    np.savez_compressed('nematic.npz', time=t.time/1000, S=S)
    plt.figure()
    plt.plot(t.time/1000, S, linewidth=2)
    plt.ylabel('Nematic Order Parameter', fontsize=14)
    plt.xlabel('Time (ns)', fontsize=14)
    plt.ylim(-0.1, 1)

    equil = pymbar.timeseries.detectEquilibration(S)[0]  # frame at which S is considered to be equilibrated

    equilibrated_frames = t.n_frames - equil
    print('Equilibration detected after %s ns' % (t.time[equil]/1000))

    if not args.blocks:
        tau = int(pymbar.timeseries.integratedAutocorrelationTime(S[equil:]))  # find the correlation time for equilibrated data
        if tau == 0:
            tau = 1
        nblocks = int(S[equil:].shape[0] / tau)
    else:
        nblocks = args.blocks
        tau = int(S[equil:].shape[0] / nblocks)

    print('Trajectory divided into %s blocks of length %s ns' % (int(nblocks), tau))

    S_avg = np.mean(S[equil:])

    data = S[-nblocks*tau:]
    block_averages = np.zeros([nblocks])
    for b in range(nblocks):
        block_averages[b] = np.mean(data[b*tau:(b + 1)*tau])
    S_std = np.std(block_averages)

    print('Nematic order parameter: %0.3f +/- %0.3f' % (S_avg, S_std))

    h = np.zeros([nT - equil, args.bins])
    for i in range(equil, nT):
        h[i - equil, :], bin_edges = np.histogram(a[i, :], bins=args.bins, range=(0, 1))

    # Boostrap to get uncertainties in histogram bars
    hboot = np.zeros([args.nboot, nblocks, args.bins])  # average histogram from bootstrap trial
    for b in range(args.nboot):
        boot = np.zeros([nblocks, args.bins])
        for n in range(nblocks):
            ndx = np.random.randint(n*tau, (n+1)*tau, tau)
            trials = h[ndx, :]  # histogram from random time points within block
            hboot[b, n, :] = np.mean(trials, axis=0)

    meanboot = np.zeros([nblocks, args.bins])
    for n in range(nblocks):
        meanboot[n, :] = np.mean(hboot[:, n, :], axis=0)

    bin_centers = [bin_edges[i] - ((bin_edges[i] - bin_edges[i - 1]) / 2) for i in range(args.bins)]
    bin_width = 1 / args.bins

    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.text(0.4, 0.7, 'Nematic Order Parameter:\n %0.3f $\pm$ %0.3f' % (S_avg, S_std), horizontalalignment='center',
             verticalalignment='center', transform=ax.transAxes, fontsize=18)
    plt.bar(bin_centers, np.mean(h, axis=0)/((a.shape[1])*bin_width), bin_width,
            yerr=np.std(meanboot, axis=0)/((a.shape[1]*bin_width)))
    # plt.bar(bin_centers, np.mean(h, axis=0)/((a.shape[1])*bin_width), bin_width,
    #         yerr=std_boot/((a.shape[1]*bin_width)))
    # plt.bar(bin_centers, np.mean(h[equil:, :], axis=0), bin_width)
    ax.tick_params(axis='both', labelsize=14)
    plt.xlabel('$cos(\Theta)$', fontsize=14)
    plt.ylabel('Frequency (%)', fontsize=14)
    plt.tight_layout()
    plt.savefig('nematic_order.png')
    plt.show()