#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import os
import tqdm
import matplotlib.pyplot as plt
from LLC_Membranes.analysis import hbonds

script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='wiggle.trr', type=str, help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate file')
    parser.add_argument('-p', '--top', default='topol.top', type=str, help='Gromacs topology file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='End frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Start frame')
    parser.add_argument('-x', '--exclude_water', action='store_true', help='Exclude water while searching for hbonds')
    parser.add_argument('-r', '--residues', nargs='+', default=['HII'], help='Residues to include in h-bond search')
    parser.add_argument('-a', '--atoms', action='append', nargs='+', help='Atoms for to include for each '
                        'residue. Each list of atoms must be passed with a separate -a flag for each residue in'
                        'args.residues')
    parser.add_argument('-d', '--distance', default=.3, help='Maximum distance between acceptor and donor atoms')
    parser.add_argument('-angle', '--angle_cut', default=20, help='Maximum DHA angle to be considered an H-bond')
    parser.add_argument('-nmon', default=5, type=int, help='Number of monomers per layer')

    args = parser.parse_args()

    return args


def different_layers(x, nmon):
    """
    Test if monomers a pair of monomers are in separate layers (based on initial configuration build procedure)
    :param x : length 2 list of monomer numbers (calculated using hbonds.System.number_residues)
    :param nmon: number of monomers per layer (int)
    :return: True / False (bool)
    """

    different = False
    layer = [a // nmon for a in x]
    if layer[0] != layer[1]:
        different = True

    return different


def running_average(mylist, N):

    cumsum, moving_aves = [0], []

    for i, x in enumerate(mylist, 1):
        cumsum.append(cumsum[i-1] + x)
        if i >= N:
            moving_ave = (cumsum[i] - cumsum[i-N])/N
            moving_aves.append(moving_ave)

    return moving_aves


if __name__ == "__main__":

    args = initialize()

    # workaround for argparse. If default value is set, it is always included in the list with action='append'
    if not args.atoms:
        args.atoms = ['O3', 'O4']  # a default value

    sys = hbonds.System(args.traj, args.gro, args.top, begin=args.begin, end=args.end, exclude_water=args.exclude_water)

    for i, r in enumerate(args.residues):
        sys.set_eligible(r, args.atoms[i])

    water_numbers = sys.number_water_molecules()
    res_numbers, nres = sys.number_residues('HII')

    sys.identify_hbonds(args.distance, args.angle_cut)

    # single = np.zeros(sys.t.n_frames)
    nlayers = 20
    pores = 4
    intra = np.zeros([pores*nlayers, sys.t.n_frames])  # h-bonds within layer
    inter = np.zeros([pores*nlayers, sys.t.n_frames])  # h-bonds between layers
    # above_and_below = np.zeros(sys.t.n_frames)
    for i in range(sys.t.n_frames):
        # check to see which hydrogen share bonds between monomer and if those monomers are in different layers
        n = 0
        donors = list(sys.hbonds[i][0, :])  # [D, H, A, angle]
        double_donors = list(set(int(x) for x in donors if donors.count(x) > 1))
        donors = np.array(donors)
        interlayer = []
        intralayer = []
        for d in double_donors:
            ndx = np.where(donors == d)[0]
            x = [res_numbers[int(x)] for x in sys.hbonds[i][2, ndx] if int(x) in res_numbers]
            if len(x) == 2:
                layer = min([a // args.nmon for a in x])
                if different_layers(x, args.nmon):
                    inter[layer, i] += 1
                    # n += 1
                    interlayer.append(x)
                else:
                    intra[layer, i] += 1

        # check if monomers sharing hbonds across layers, share another hbond with a third layer
        # N = 0
        # for l in range(len(interlayer)):
        #     for j in range(len(interlayer)):
        #         if interlayer[l][0] in interlayer[j] and j != l:
        #             ndx = interlayer[j].index(interlayer[l][0])
        #             # check that monomer isn't sharing hbond with 2 monomers in the same layer
        #             if different_layers([interlayer[l][1], interlayer[j][1 - ndx]], args.nmon):
        #                 N += 1
        #         if interlayer[l][1] in interlayer[j] and j != l:
        #             ndx = interlayer[j].index(interlayer[l][1])
        #             if different_layers([interlayer[l][0], interlayer[j][1 - ndx]], args.nmon):
        #                 N += 1
        # N /= 2  # avoid double counting
        # single[i] = n
        # above_and_below[i] = N

    fig1, ax1 = plt.subplots(2, 2)
    for i in range(2):
        for j in range(2):
            ax1[i, j].plot(np.linspace(1, nlayers, nlayers), np.mean(inter[(2*i + j)*nlayers:(2*i + j + 1)*nlayers, :], axis=1), linewidth=2)
            ax1[i, j].text(.5, .85, 'Pore %d' % (1 + j + 2 * i), horizontalalignment='center', transform=ax1[i, j].transAxes, fontsize=14)
            ax1[i, j].set_ylim(0, 1)

    ax1[1, 0].set_xlabel('Layer', fontsize=14)
    ax1[1, 1].set_xlabel('Layer', fontsize=14)
    ax1[0, 0].set_ylabel('Average shared\n hbonds', fontsize=14)
    ax1[1, 0].set_ylabel('Average shared\n hbonds', fontsize=14)

    plt.tight_layout()
    plt.savefig('pore_hbonds.png')

    # pad between data points with zeros
    block_length = 2  # insert this many zeros after each entry in layers
    pores = 4
    layers = np.linspace(0, pores*nlayers - 1, pores*nlayers*(block_length + 1) - block_length)
    binsize = (layers[1] - layers[0])

    ft = np.zeros([pores, int(nlayers*(block_length + 1) - block_length), sys.t.n_frames])

    for j in range(inter.shape[1]):
        for k in range(4):

            avg = np.copy(inter[k*nlayers:(k+1)*nlayers, j])

            for i in range(avg.size - 1):
                avg = np.insert(avg, 1 + i + i*block_length, np.zeros(block_length))

            avg -= avg.mean()

            fft = np.abs(np.fft.fft(avg)) ** 2

            # filter out certain lower frequencies
            # W = np.fft.fftfreq(avg.size, d=binsize)
            # cut_f_signal = np.copy(fft)
            # cut_f_signal[(W >= 0.9)] = 0
            # cut_f_signal[(W <= -0.9)] = 0
            # inv = np.fft.ifft(cut_f_signal)
            # fft = np.abs(np.fft.fft(inv)) ** 2
            # #####################################

            freq_x = np.fft.fftfreq(avg.size, d=binsize)
            ndx = np.argsort(freq_x)
            freq_x = freq_x[ndx]
            ft[k, :, j] += fft[ndx]

    # plt.plot(np.linspace(0, pores*nlayers - 1, pores*nlayers), np.mean(inter, axis=1))
    fig, ax = plt.subplots(2, 2)
    for i in range(2):
        for j in range(2):
            ax[i, j].plot(freq_x, np.mean(ft[2*i + j, :, :], axis=1), linewidth=2)
            ax[i, j].text(.5, .85, 'Pore %d' % (1 + j + 2 * i), horizontalalignment='center', transform=ax[i, j].transAxes, fontsize=14)
            ax[i, j].set_ylim(0, 7)
            ax[i, j].set_xlim(-1.2, 1.2)

    ax[1, 0].set_xlabel('Frequency (layers$^{-1}$)', fontsize=14)
    ax[1, 1].set_xlabel('Frequency (layers$^{-1}$)', fontsize=14)
    ax[0, 0].set_ylabel('Intensity', fontsize=14)
    ax[1, 0].set_ylabel('Intensity', fontsize=14)

    plt.tight_layout()
    plt.savefig('pore_hbonds_ft.png')

    # plt.xlabel('Layer')
    # plt.ylabel('Instensity')
    # plt.figure()
    # plt.plot(freq_x, np.mean(ft, axis=1))
    plt.show()
    exit()

    # these both work out to be equal to the number of monomer residues
    total_possible_interlayer = nres  # number of possible pairs that share a water molecule
    total_possible_combos = nres  # number of possible triplets with a continuous hbond connection

    p_interlayer = single / total_possible_interlayer
    p_combos = above_and_below / total_possible_combos

    p = (p_interlayer**2 - p_combos) / p_interlayer

    plt.plot(sys.t.time, p)

    plt.figure()
    plt.bar(sys.t.time, p_combos, sys.t.time[1] - sys.t.time[0])

    plt.figure()
    plt.bar(sys.t.time, p_interlayer, sys.t.time[1] - sys.t.time[0])
    N = 50
    plt.plot(sys.t.time[(N-1):], running_average(p_interlayer, N), linewidth=2)

    plt.show()
