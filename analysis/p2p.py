#!/usr/bin/env python


from __future__ import division
from __future__ import print_function
import numpy as np
import math
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import griddata
from matplotlib import animation
import argparse
from pymbar import timeseries
import random as ran
import mdtraj as md
from scipy.optimize import curve_fit
from scipy import spatial
import tqdm
from scipy.misc import factorial


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the pore-to-pore distance from an MD trajectory of an '
                                                 'HII LLC membrane')

    # Trajectory loading and slicing
    parser.add_argument('-t', '--input', default='wiggle.trr', type=str, help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate configuration file (.gro, .pdb)')
    parser.add_argument('--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('--end', default=-1, type=int, help='Frame to stop calculations')
    parser.add_argument('--skip', default=1, type=int, help='Sample every nth frame')

    # System-dependent parameters
    parser.add_argument('-p', '--pores', default=4, help='Number of pores in unit cell')
    parser.add_argument('-c', '--component', default='NA', nargs='+', help='Name of monomer component used to track'
                        'pore positions. This script will calculate the average coordinate of each in each pore. '
                        'Special predefined cases include: "tails", "Tails", "benzene", "Benzene", "Head Groups, '
                        '"tail_ends", "tail_fronts", "Sodium". These are specific to Na-GA3C11')
    parser.add_argument('-E', '--equil', default='auto', help='Frame number where system is equilibrated. "auto" will '
                        'use pymbar.timeseries.DetectEquilibration to determine which frame to start at. It is often '
                        'worth double checking its choice manually')
    parser.add_argument('-x', '--exclude', default=[4], nargs='+', help='Index of p2p distance to exclude as written in'
                        'the list: ["1-2", "1-3", "1-4", "2-3", "2-4", "3-4"] ')
    parser.add_argument('--auto_exclude', action="store_true", help="Specifying this will override args.exclude and "
                        "decide which pore-to-pore distance to exclude automatically by dropping the highest value")
    parser.add_argument('-b', '--nboot', default=2000, help='Number of bootstrap trials for generating statistics')

    # plotting details
    parser.add_argument('--plot_every', default=1, type=int, help='Plot every n frames')
    parser.add_argument('--noshow', help='Specify this flag to prevent the plot from showing', action="store_true")
    parser.add_argument('--save', action="store_true", help='Save the output plot')
    parser.add_argument('--plot_avg', action="store_true", help='Plot average p2p distance at each frame')
    parser.add_argument('--plot_std', action="store_true", help='Plot average p2p distance at each frame')

    # special cases
    parser.add_argument('-T', help='Use for concatenated trajectories of simulations run at multiple '
                        'temperatures. Supply the name of index file (ascii format) with temperatures and times where '
                        'the switch occurs. On each line supply the info with the format temp:time')
    parser.add_argument('--buffer', default=0, type=float, help='Fraction (of membrane thickness) of top and bottom of'
                        'membrane to exclude from p2p calculations. Useful if you solvate both sides with water and'
                        'ions float into solution.')

    return parser


def restrict_atoms(t, component):
    """ Restrict trajectory of coordinates to selected atoms

    :param t: mdtraj trajectory object
    :param component: names of components to be used for tracking pore centers.

    :type t: mdtraj.core.trajectory.Trajectory
    :type component: str or list
    :return:
    """

    # Check for certain special arguments
    if component == 'tails' or component == 'Tails':
        atoms = ['C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21',
                 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35',
                 'C36', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48']
    elif component == 'benzene' or component == 'Benzene':
        atoms = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
    elif component == 'Head Groups':
        atoms = ['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O', 'O1', 'O2', 'O3', 'O4']
    elif component == 'tail_ends':
        atoms = ['C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C24', 'C25', 'C26', 'C27',
                 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45',
                 'C46', 'C47', 'C48']
    elif component == 'tail_fronts':
        atoms = ['C7', 'C8', 'C9', 'C21', 'C22', 'C23', 'C35', 'C36', 'C37']
    elif component == 'Sodium':
        atoms = ['NA']
    else:
        atoms = component

    if component != 'sys':

        atoms_to_keep = [a.index for a in t.topology.atoms if a.name in atoms and a.residue.name != 'HOH']
        pos = t.xyz[:, atoms_to_keep, :]

    else:
        # NOTE: if you use 'sys', use gmx trjconv -f *.trr -s *.tpr -pbc atom
        #                        then gmx trjconv -f *.trr -s *.tpr -pbc whole
        # This will make sure everything is in the box and whole.
        pos = t.xyz

    return pos


def avg_pore_loc(npores, pos, buffer=0):
    """ Calculate average pore location for each pore at each frame

    :param no_pores: the number of pores in the unit cell
    :param pos: the coordinates of the component(s) which you are using to locate the pore centers
    :param buffer: fraction (of membrane thickness) of top and bottom of membrane to exclude from p2p calculations

    :type no_pores: int
    :type pos: numpy.ndarray, shape(ncomponents, 3) or numpy.ndarray, shape(nframes, ncomponents, 3)
    :type buffer: float

    :return: numpy array containing the x, y coordinates of the center of each pore at each frame
    """

    # Find the average location of the pores w.r.t. x and y

    if len(pos.shape) == 3:  # multiple frames

        nT = np.shape(pos)[0]
        comp_ppore = np.shape(pos)[1] // npores

        p_center = np.zeros([2, npores, nT])

        for i in range(nT):
            zmax = np.amax(pos[i, :, 2])  # maximum z value for this frame
            zmin = np.amin(pos[i, :, 2])  # minimum z value for this frame
            thick = zmax - zmin
            zmax -= buffer*thick
            zmin += buffer*thick
            for j in range(npores):
                count = 0
                for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                    if zmax >= pos[i, k, 2] >= zmin:
                        p_center[:, j, i] += pos[i, k, :2]
                        count += 1
                p_center[:, j, i] /= count  # take the average

    elif len(pos.shape) == 2:  # single frame

        comp_ppore = pos.shape[0] // npores

        p_center = np.zeros([2, npores])

        for j in range(npores):
            for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                p_center[:, j] += pos[k, :2]
            p_center[:, j] /= comp_ppore

    else:
        return 'Please use a position array with valid dimensions'
        exit()

    return p_center


def p2p(p_centers, distances):
    """ Calculate all pairwise pore-to-pore distances (4 pores is the only number of pores implemented currently)

    :param p_centers: the x, y locations of the pore centers
    :param distances: the number of distinct distances between pores

    :type p_centers: numpy.ndarray, shape(p_centers.shape[0], p_centers.shape[1], 2)
    :type distances: int

    :return: (np.ndarray, shape(p_centers.shape[0], distances, 2) All frame-by-frame pore-to-pore distances
    """

    nT = np.shape(p_centers)[2]
    p2ps = np.zeros([distances, nT])  # distances in the order 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
    for i in range(nT):
        # This could potentially be improved by figuring out an empirical formula for number of p2p distances of
        # interest as a function of pores. Then only return the bottom x number of p2p distances.
        p2ps[0, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 1, i])
        p2ps[1, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 2, i])
        p2ps[2, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 3, i])
        p2ps[3, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 2, i])
        p2ps[4, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 3, i])
        p2ps[5, i] = np.linalg.norm(p_centers[:, 2, i] - p_centers[:, 3, i])

    return p2ps


def p2p_stats(p2ps, exclude, nboot, equil):
    """ Calculate the average and spread of pore-to-pore distances

    :param p2ps: all of the pore-to-pore distances
    :param exclude: exclude a certain pore-to-pore interaction (such as the distance from pores on opposite sides of a
           parallelogram as is the case in a hexagonal system). If 'auto', then exclude will be automatically selected
    :param nboot: number of bootstrap trials to use when generating statistics
    :param equil: the trajectory frame at which to start generating statistics. Care about this parameter if you
           choose to detect equilibration manually. Otherwise 'auto' will use pymbar to find it for you

    :type p2ps: numpy.ndarray, shape(nframes, np2p_distances)
    :type exclude: int or str
    :type nboot: int
    :type equil: int or str

    :return: the average and standard deviation of pore to pore distances
    """

    nboot = int(nboot)
    nT = p2ps.shape[1]
    ndist = p2ps.shape[0]
    ndist -= len(exclude)

    # Get rid of the excluded trajectory
    p2p_new = np.zeros([ndist, nT])
    count = 0
    for i in range(ndist + len(exclude)):
        if i not in exclude:
            p2p_new[count, :] = p2ps[i, :]
            count += 1

    p2ps = p2p_new

    # Find the frame at which the system is equilibrated
    if equil == 'auto':
        ts = []
        for pore in range(ndist):
            ts.append(timeseries.detectEquilibration(p2ps[pore, :])[0])
        t = int(max(ts))  # use the max equil frame to ensure all pores are equilibrated
    else:
        t = int(equil)

    # Find the autocorrelation time for each pore - i.e. the time it takes for samples to become uncorrelated
    taus = []
    for i in range(ndist - len(exclude)):
        tau = timeseries.integratedAutocorrelationTime(p2ps[i, t:])
        taus.append(tau)

    print('Maximum Autocorrelation Time: %s frames' % max(taus))
    tau = int(np.ceil(max(taus)))  # use the max again to ensure all trajectories are independent. np.ceil e

    if tau == 0:
        tau = 1

    ind_trajectories = (nT - t) // tau  # the number of independent trajectories
    print('%s Independent Trajectories' % ind_trajectories)

    trajectories = np.zeros([ndist, ind_trajectories, tau])  # Create a new array to hold all the trajectories

    # fill up trajectory array
    for i in range(ndist):
        for j in range(ind_trajectories):
            # trajectories[i*ind_trajectories + j, :] = p2ps[i, (t + j*tau):(t + (j + 1)*tau)]
            trajectories[i, j, :] = p2ps[i, (t + j*tau):(t + (j + 1)*tau)]

    # bootstrap to get statistics
    p2p_boot = np.zeros([nboot, ndist, ind_trajectories*tau])  # a bunch of full trajectories assembled from independent trajectories
    avg_trials = np.zeros([nboot, ndist])
    for i in range(nboot):
        for k in range(ind_trajectories):
            T = ran.randrange(0, ind_trajectories)  # pick a random trajectory from all the independent trajectories
            p2p_boot[i, :, k*tau:(k+1)*tau] = trajectories[:, T, :]
        avg_trials[i, :] = np.mean(p2p_boot[i, :, :], axis=1)  # Average value of each pore for each trial

    average_distances = np.mean(avg_trials, axis=0)
    avg = np.mean(average_distances)

    var = 0
    for i in average_distances:
        var += (i - avg)**2

    var /= 4  # number of independent distances
    std = var ** 0.5  # take square root to get std from variance

    return avg, std, t


def parse_txt(txt, points, times, std=False):
    """ Find out where to place lines in plot

    :param txt: a text file with pertinent information. Temperatures which the simulation was run at and at what time
    the simulation is at that temperature should be of the form Temperature:time (i.e. at `time` ps the simulation
    switches to `Temperature`). That section should be donoted by a header of the form [ T ]. You can also include a
    benchmark value which will be plotted as a straight line with dotted lines representing the error bounds. Start that
    section with the header [ bench ] with subsequent lines of the form value+/-error. Sections should be
    separated by a blank line
    :param points:

    :type txt: str

    :return: vertical lines where temperature switches occur (with T labels) and horizontal lines for the benchmark
    value with dotted horizontal lines representing the error bounds
    """

    with open(txt, 'r') as f:

        a = []
        for line in f:
            a.append(line)

    T = 0
    while a[T].count('[ T ]') == 0:
        T += 1
    T += 1
    temp = []
    time = []
    while a[T] != '\n':
        info = a[T].split(':')
        temp.append(int(info[0]))
        time.append(float(info[1]))
        T += 1

    for i, T in enumerate(temp):
        x = np.array([time[i], time[i]])
        y = np.array([0, 5])
        plt.plot(x, y, '--', linewidth=3, color='r', scaley=False)
        plt.annotate('T = %s' % T, xy=(time[i], np.mean(points)), xytext=(time[i] + 10, np.amax(points)))

    if std:
        std_equil = np.zeros([4, len(time) - 1])
        for i in range(len(time) - 1):
            begin = np.where(times == time[i])[0][0]
            end = np.where(times == time[i + 1])[0][0]
            slice = p2ps[:, begin:end]
            p2p_avg, p2p_std, equil = p2p_stats(slice, exclude, '%s' % args.nboot, '%s' % args.equil)
            std_equil[:, i] = [p2p_avg, p2p_std, equil + begin, end]

        for i in range(std_equil.shape[1]):
            avg = std_equil[0, i]
            std = std_equil[1, i]
            x = np.array([times[int(std_equil[2, i])], times[int(std_equil[3, i])]])
            upperbound = np.array([avg + std, avg + std])
            lowerbound = np.array([avg - std, avg - std])
            plt.fill_between(x, upperbound, lowerbound, alpha=0.5, color='orange')
            # plt.plot(x, upperbound, '--', linewidth=1, color='orange')
            # plt.plot(x, lowerbound, '--', linewidth=1, color='orange')

    B = 0
    while a[B].count('[ bench ]') == 0:
        B += 1
    B += 1

    benchval = []
    bencherr = []
    y_low = [np.amin(points)]
    y_high = [np.amax(points)]

    while a[B] != '\n' and B < len(a):
        info = a[B].split('+/-')
        benchval.append(float(info[0]))
        bencherr.append(float(info[1]))
        B += 1

    for i, v in enumerate(benchval):
        if i < len(temp) - 1:
            x = np.array([time[i], time[i + 1]])
        else:
            x = np.array([time[i], (time[i] + 1)*1000000])
        y = np.array([v, v])
        y_low.append(v - bencherr[i])
        y_high.append(v + bencherr[i])
        upperbound = np.array([v + bencherr[i], v + bencherr[i]])
        lowerbound = np.array([v - bencherr[i], v - bencherr[i]])
        plt.fill_between(x, upperbound, lowerbound, alpha=0.5, color='g')
        plt.plot(x, y, linewidth=2, color='k', scalex=False)
        # plt.plot(x, upperbound, '--', linewidth=1, color='g')
        # plt.plot(x, lowerbound, '--', linewidth=1, color='g')

    ybounds = np.array([min(y_low)*.99, max(y_high)*1.01])
    xbounds = np.array([-3, max(times)])

    return xbounds, ybounds


if __name__ == '__main__':

    args = initialize().parse_args()  # parse the args

    t = md.load(args.input, top=args.gro)[args.begin:args.end:args.skip]

    pos = restrict_atoms(t, args.component)  # convenience function
    nT = np.shape(pos)[0]

    tot_atoms = np.shape(pos)[1]
    n_pores = int(args.pores)  # number of pores
    comp_ppore = tot_atoms // n_pores

    p_centers = avg_pore_loc(n_pores, pos, args.buffer)

    distances = 6  # number of p2p distances to calculate. The following algorithm isn't smart enough for >6 yet
    p2ps = p2p(p_centers, distances)

    if args.auto_exclude:
        means = np.zeros([p2ps.shape[0]])
        for i in range(p2ps.shape[0]):
            means[i] = np.mean(p2ps[i, :])
        exclude = [np.argmax(means)]
    else:
        exclude = [int(i) for i in args.exclude]

    p2p_avg, p2p_std, equil = p2p_stats(p2ps, exclude, '%s' % args.nboot, '%s' % args.equil)
    print('Equilibration detected after %d ns' % (t.time[equil] / 1000))
    print('Average Pore to Pore distance: %.3f' % p2p_avg)
    print('Standard Deviation of Pore to Pore distances: %.3f' % p2p_std)

    labels = ['1-2', '1-3', '1-4', '2-3', '2-4', '3-4']
    labels = ['1-4', '1-0', '0-2', '0-4', '0-3', '2-3', '2-5', '5-3', '3-4', '5-6', '3-7', '6-7', '4-7', '7-8',
              '4-8', '3-6']

    if args.plot_avg:

        fig = plt.figure(1)
        ax1 = fig.add_subplot(111)
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Distance between pores (nm)', fontsize=14)
        ax1.tick_params(axis='both', labelsize=14)

        avg = np.zeros(p2ps[0, :].shape)
        n = 0
        for i in range(distances):
            if i not in exclude:
                avg += p2ps[i, :]
                n += 1

        avg /= n

        ax1.plot(t.time[::args.plot_every]/1000, avg[::args.plot_every], linewidth=2)

        if args.T:

            xbounds, ybounds = parse_txt(args.T, avg[::args.plot_every], t.time[::args.plot_every]/1000, std=args.plot_std)
            plt.ylim(ybounds)
            plt.xlim(xbounds)

        plt.tight_layout()

    else:
        fig = plt.figure(1)
        ax1 = fig.add_subplot(111)

        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Distance between pores (nm)', fontsize=14)
        ax1.tick_params(axis='both', labelsize=14)

        for i in range(distances):
            if i not in exclude:
                plt.plot(t.time[::args.plot_every], p2ps[i, ::args.plot_every], label='%s' % labels[i])

        plt.tight_layout()

    if args.save:
        plt.savefig('p2p.png')
    if not args.noshow:
        plt.show()
