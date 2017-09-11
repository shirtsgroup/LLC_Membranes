#! /usr/bin/env python

import argparse
import subprocess
import mdtraj as md
import Structure_char
import numpy as np
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Using three output files from an MD simulation, run all analysis scripts')

    parser.add_argument('-t', '--traj', default='wiggle.trr', type=str, help='Name of trajectory file (.xtc or .trr)')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-tpr', '--tpr', default='wiggle.tpr', type=str, help='Name of atomic level input file (.tpr)')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Start frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop calculations')
    parser.add_argument('--skip', default=1, type=int, help='Usage: --skip n . Sample every nth frame')
    parser.add_argument('-a', '--atoms', nargs='+', type=str, default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'],
                        help='Name of atoms to calculate z distribution function with respect to')
    parser.add_argument('--noxrd', action="store_true", help='Do not simulate x-ray diffraction pattern')
    parser.add_argument('-xf', '--xrd_frames', default=-50, help='Start frame for xrd simulation')
    parser.add_argument('--plot_pores', action="store_true", help='Plot trajectories of individual pore-to-pore distances')

    args = parser.parse_args()

    return args

# def tex_figure():


if __name__ == "__main__":

    args = initialize()

    # Modify the trajectory so all monomers stay whole
    # echo = subprocess.Popen(["echo", "0"], stdout=subprocess.PIPE)
    # ps = subprocess.Popen(["gmx", "trjconv", "-f", "%s" % args.traj, "-o", "traj_whole.xtc", "-s", "%s" % args.tpr,
    #                       "-pbc", "whole"], stdin=echo.stdout)
    # ps.wait()
    t = md.load(args.traj, top=args.gro)[args.begin:args.end:args.skip]
    nT = t.n_frames

    if not args.noxrd:

        subprocess.call(["main_gromacs.py", "-top", "%s" % args.gro, "-traj", "traj_whole.xtc", "--lcscale", "1.2",
                         "-fi", "%s" % args.xrd_frames])

    # Pore-to-pore distances

    # constants
    n_pores = 4  # number of pores
    distances = 6  # number of p2p distances to calculate. My algorithm isn't smart enough for anything but six yet

    pos = Structure_char.restrict_atoms(t, 'NA')

    p_centers = Structure_char.avg_pore_loc(n_pores, pos, 0)

    p2ps = Structure_char.p2p(p_centers, 6)

    means = np.zeros([p2ps.shape[0]])
    for i in range(p2ps.shape[0]):
        means[i] = np.mean(p2ps[i, :])
    exclude = [np.argmax(means)]

    p2p_avg, p2p_std, equil = Structure_char.p2p_stats(p2ps, exclude, 2000, 'auto')
    print 'Equilibration detected after %d ns' % (t.time[equil] / 1000)
    print 'Average Pore to Pore distance: %.3f' % p2p_avg
    print 'Standard Deviation of Pore to Pore distances: %.3f' % p2p_std

    labels = ['1-2', '1-3', '1-4', '2-3', '2-4', '3-4']

    if not args.plot_pores:
        plt.figure(1)
        avg = np.zeros(p2ps[0, :].shape)
        n = 0
        for i in range(distances):
            if i not in exclude:
                avg += p2ps[i, :]
                n += 1

        avg /= n

        plt.plot(t.time, avg)

    else:
        plt.figure(1)
        for i in range(distances):
            if i not in exclude:
                plt.plot(t.time, p2ps[i, :], label='%s' % labels[i])

    plt.title('Pore to Pore Distance Equilibration')
    plt.ylabel('Distance between pores (nm)')
    plt.xlabel('Time (ps)')
    if args.plot_pores:
        plt.legend(loc=1, fontsize=18)

    plt.savefig('p2p.png')

