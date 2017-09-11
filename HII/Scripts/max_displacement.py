#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
import tqdm


def initialize():

    parser = argparse.ArgumentParser(description='Calculate maximum z displacement of an atom or molecule')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Name of trajectory file (.xtc or .trr)')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop calculations')
    parser.add_argument('--skip', default=1, type=int, help='Usage: --skip n . Sample every nth frame')
    parser.add_argument('-a', '--atoms', nargs='+', type=str, default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'],
                        help='Name of atoms to track')
    parser.add_argument('--time', action="store_true", help='Plot the average max displacement against time')
    parser.add_argument('--save', action="store_true", help='Save output figures')

    args = parser.parse_args()

    return args


def time_displacement(pos):
    """
    :param pos: z position of each atom versus time [nT, natoms]
    :return: average max displacement versus time
    """

    nT = pos.shape[0]
    natoms = pos.shape[0]

    td = np.zeros([nT])
    itau = 1

    for itau in tqdm.tqdm(range(1, nT)):
        ncount = natoms*(nT - itau)
        for n in range(natoms):
            xn = pos[:, n]
            for t in range(nT - itau):
                xo = np.max(xn[t:(t + itau)]) - np.min(xn[t:(t + itau)])
                td[itau] += xo

        td[itau] /= ncount
        itau += 1

    return td


# def time_displacement(pos):
#     """
#     :param pos: z position of each atom versus time [nT, natoms]
#     :return: average max displacement versus time
#     """
#
#     nT = pos.shape[0]
#     natoms = pos.shape[0]
#
#     tdisp = np.zeros([nT, natoms])
#
#     for i in range(natoms):
#         for t in range(nT):
#             tdisp[t, i] = np.max(pos[:(t+1), i]) - np.min(pos[:(t+1), i])
#
#     tavg = np.zeros([nT])
#     for i in range(nT):
#         tavg[i] = np.mean(tdisp[i, :])
#
#     return tavg


def displacements(pos):
    """
    :param pos: z position of each atom versus time [nT, natoms]
    :return: max displacement for each atom
    """

    natoms = pos.shape[1]
    disp = np.zeros([natoms])

    for a in range(natoms):
        disp[a] = np.max(pos[:, a]) - np.min(pos[:, a])

    return disp

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[args.begin:args.end:args.skip]
    nT = t.n_frames

    box = t.unitcell_vectors
    L = np.min(box[:, 2, 2])  # average z length of unit cell

    keep = [a.index for a in t.topology.atoms if a.name in args.atoms]

    pos = t.atom_slice(keep).xyz[:, :, 2]  # z positions only

    if args.time:

        d = time_displacement(pos)

        plt.plot(t.time, d)
        plt.xlabel('Time (ps)')
        plt.ylabel('Average maximum displacement')
        if args.save:
            plt.savefig('average_max_displacement.png')

    else:

        d = displacements(pos)

        d = np.ma.masked_greater(d, L)

        a = np.argmax(d)
        mon = np.floor(a/len(args.atoms))
        print(mon)
        print('Watch atoms %s to %s' %(int(mon*137), int((mon + 1)*137)))
        print(d[a])

        plt.hist(d.compressed(), bins=25)
        plt.xlabel('Maximum displacement (nm)')
        plt.ylabel('Count')
        if args.save:
            plt.savefig('max_displacement.png')

    plt.show()