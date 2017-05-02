#!/usr/bin/env python

# This is the revamped version of Cylindricity.py

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


def initialize():
    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--input', default='wiggle.trr', help = 'Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help = 'Some kind of configuration file that mdtraj needs'
                                                                    'Can be a .pdb as well')
    parser.add_argument('-n', '--no_monomers', default=6, help = 'Number of Monomers per layer')
    parser.add_argument('-a', '--atoms', default=137, help = 'Number of atoms per monomer')
    parser.add_argument('-l', '--layers', default=20, help = 'Number of layers in each pore')
    parser.add_argument('-p', '--pores', default=4, help = 'Number of Pores')
    parser.add_argument('-c', '--component', default='NA', help = 'Counterion used to track pore positions')
    parser.add_argument('-f', '--start_frame', default=0, type=int, help = 'Frame number to start reading trajectory at')
    parser.add_argument('-s', '--layer_distribution', default='uniform', help = 'The distribution of monomers per layer')
    parser.add_argument('-L', '--alt_1', default=6, help = 'Monomers per layer for the first type of alternating layer')
    parser.add_argument('-A', '--alt_2', default=8, help = 'Monomers per layer for the second type of alternating layer')
    parser.add_argument('-d', '--direction', default='z', help = 'Axis along which to measure a component density')
    parser.add_argument('-S', '--slices', default='250', help = 'Number of slices to descretize chosen axis direction into')
    parser.add_argument('-C', '--cmap', default='Blues', help = 'Color Scheme for heat map')
    parser.add_argument('-m', '--lc', default='HII', help = 'Type of liquid crystal monomer being used')
    parser.add_argument('-E', '--equil', default='auto', help = 'Frame number where system is equilibrated. "auto" will '
                        'use pymbar.timeseries.DetectEquilibration to determine which frame to start at. It is worth '
                        'double checking its choice manually')
    parser.add_argument('-x', '--exclude', default=[4], nargs='+', help = 'Which pore-to-pore distance to exclude - pass the index of'
                                                            'the pore-to-pore distance as written in the list: '
                                                            '["1-2", "1-3", "1-4", "2-3", "2-4", "3-4"] ')
    parser.add_argument('-b', '--nboot', default=2000, help = 'Number of bootstrap trials')
    parser.add_argument('--noshow', help='Specify this flag to prevent the plot from showing', action="store_true")
    parser.add_argument('--auto_exclude', action="store_true", help="Specifying this will decide which pore-to-pore"
                                                    "distance to exclude automatically by dropping the highest value")
    parser.add_argument('--save', action="store_true", help='Save the output plot')

    args = parser.parse_args()
    return args


def avg_pore_loc(npores, pos, natoms):
    """
    :param no_pores: the number of pores in the unit cell
    :param pos: the coordinates of the component(s) which you are using to locate the pore centers
                      (numpy array with dimensions: [no frames, no components, xyz coordinates, ] or just
                      [no components, xyz coordinates] for a single frame)
    :return: numpy array containing the x, y coordinates of the center of each pore at each frame
    """

    # Find the average location of the pores w.r.t. x and y

    if len(pos.shape) == 3:  # multiple frames

        nT = np.shape(pos)[0]
        comp_ppore = np.shape(pos)[1] / npores

        p_center = np.zeros([2, npores, nT])

        for i in range(nT):
            for j in range(npores):
                for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                    p_center[:, j, i] += pos[i, k, :2]
                p_center[:, j, i] /= comp_ppore  # take the average

    elif len(pos.shape) == 2:  # single frame

        comp_ppore = pos.shape[1] / npores
        p_center = np.zeros([2, npores])

        for j in range(npores):
            for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                p_center[:, j] += pos[:2, k]
            p_center[:, j] /= comp_ppore

    else:
        return 'Please use a position array with valid dimensions'
        exit()

    return p_center


def p2p(p_centers, distances):
    """
    :param p_centers: the x, y locations of the pore centers in the format return from avg_pore_loc()
    :param distances: the number of distinct distances between pores
    :return: all of the pore to pore distances
    """
    nT = np.shape(p_centers)[2]
    p2ps = np.zeros([distances, nT])  # distances in the order 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
    for i in range(nT):
        # So ugly ... sadness :(
        p2ps[0, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 1, i])
        p2ps[1, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 2, i])
        p2ps[2, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 3, i])
        p2ps[3, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 2, i])
        p2ps[4, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 3, i])
        p2ps[5, i] = np.linalg.norm(p_centers[:, 2, i] - p_centers[:, 3, i])
        # p2ps[0, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 4, i])
        # p2ps[1, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 1, i])
        # p2ps[2, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 2, i])
        # p2ps[3, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 4, i])
        # p2ps[4, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 3, i])
        # p2ps[5, i] = np.linalg.norm(p_centers[:, 2, i] - p_centers[:, 3, i])
        # p2ps[6, i] = np.linalg.norm(p_centers[:, 2, i] - p_centers[:, 5, i])
        # p2ps[7, i] = np.linalg.norm(p_centers[:, 5, i] - p_centers[:, 3, i])
        # p2ps[8, i] = np.linalg.norm(p_centers[:, 4, i] - p_centers[:, 3, i])
        # p2ps[9, i] = np.linalg.norm(p_centers[:, 5, i] - p_centers[:, 6, i])
        # p2ps[10, i] = np.linalg.norm(p_centers[:, 7, i] - p_centers[:, 3, i])
        # p2ps[11, i] = np.linalg.norm(p_centers[:, 6, i] - p_centers[:, 7, i])
        # p2ps[12, i] = np.linalg.norm(p_centers[:, 4, i] - p_centers[:, 7, i])
        # p2ps[13, i] = np.linalg.norm(p_centers[:, 7, i] - p_centers[:, 8, i])
        # p2ps[14, i] = np.linalg.norm(p_centers[:, 4, i] - p_centers[:, 8, i])
        # p2ps[15, i] = np.linalg.norm(p_centers[:, 3, i] - p_centers[:, 6, i])

    return p2ps


def p2p_stats(p2ps, exclude, nboot, equil):
    """
    :param p2ps: all of the pore to pore distances in the format output by p2p()
    :param exclude: int, exclude a certain pore to pore interaction (such as the distance from pores on opposite sides of a
           parallelogram as is the case in a hexagonal system)
    :param nboot: int, number of bootstrap trials to use when generating statistics
    :param equil: int, the trajectory frame at which to start generating statistics. Care about this parameter if you
           choose to detect equilibration manually. Otherwise 'auto' will use pymbar to find it for you
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

    print 'Maximum Autocorrelation Time: %s frames' % max(taus)
    tau = int(np.ceil(max(taus)))  # use the max again to ensure all trajectories are independent. np.ceil e

    if tau == 0:
        tau = 1

    ind_trajectories = (nT - t) / tau  # the number of independent trajectories

    trajectories = np.zeros([ndist, ind_trajectories, tau])  # Create a new array to hold all the trajectories

    # fill up trajectory array
    for i in range(ndist):
        for j in range(ind_trajectories):
            # trajectories[i*ind_trajectories + j, :] = p2ps[i, (t + j*tau):(t + (j + 1)*tau)]
            trajectories[i, j, :] = p2ps[i, (t + j*tau):(t + (j + 1)*tau)]

    # bootstrap to get statistics
    p2p_boot = np.zeros([ndist, nboot*ind_trajectories])  # bootstrap each pore
    for i in range(nboot):
        for j in range(ndist):
            for k in range(ind_trajectories):
                T = ran.randrange(0, ind_trajectories)  # pick a random trajectory from all the independent trajectories
                P = ran.randrange(0, tau)  # choose a random point in that trajectory
                p2p_boot[j, i*ind_trajectories + k] = trajectories[j, T, P]

    avg = np.zeros([ndist])  # find the average of each trajectory
    for i in range(ndist):
        avg[i] = np.mean(p2p_boot[i, :])

    ensemble_avg = np.mean(avg)  # the ensemble average, <x>, is the average p2p of all the trajectories

    stds = np.zeros([ndist])  # find the standard deviation for each pore : s = sqrt((<x> - x)^2) where x is the avg[i]
    for i in range(ndist):
        stds[i] = np.sqrt((ensemble_avg - avg[i])**2)

    p2p_std = np.mean(stds)  # report the average of the standard deviations

    print ' Pore to Pore Statistics calculated starting at frame {:d} ({:2.1f} percent into simulation)'.format(
            t, 100.0 * t / nT)
    return ensemble_avg, p2p_std, t


def compdensity(component, pore_centers, start, pores=4, bin_width=0.05, rmax=3.5, buffer=0.0):

    """
    :param component: the coordinates of the component(s) which you want a radial distribution of at each frame
                      (numpy array with dimensions: [3 (xyz coordinates), no components, no frames])
    :param pore_centers: a numpy array of the locations of each pore center at each trajectory frame
    :param start: the frame number at which to start calculations (should be after equilibration)
    :param step: the size of the step to take in the radial direction when measuring density (float)
    :param pores: number of pores (int) default=4
    :param bin_width: width of the bins which will show up when you plot this (float), default = 0.1 nm
    :param rmax: maximum distance from pore center to calculate density for, default = 3.5 nm
    :param buffer: percentage used to define the location of z planes between which component density will be computed,
           float, default = 0 (i.e. no buffer). Should be between 0 and 1. e.g. for 1 percent, use 0.01 as the buffer
    :return: the density of the component (comp_density) as a function the distance from the pore center (x). Also
             returns the calculated bin width for plotting
    """

    # Extract basic system information. It's important to follow the format of the component array to get it right
    n_atoms = np.shape(component)[1]  # the total number of components in a single frame
    n_ppore = tot_atoms / pores  # the total number of components in each pore
    nT = np.shape(component)[0]

    # Find the approximate max and minimum z values of the components based on the last frame
    zmax = np.max(component[-1, :, 2])
    zmin = np.min(component[-1, :, 2])
    thickness = zmax - zmin  # approximate membrane thickness

    # now find the maximum and minimum permissible z dimensions based on the buffer
    zmax_buff = zmax - thickness * buffer  # Could use buffer/2 depending on how you interpret what buffer % means
    zmin_buff = zmin + thickness * buffer

    dist_from_center = []
    # Now find the distance from the center of every atom in every frame
    for k in range(start, nT):
        for i in range(pores):
            for j in range(n_ppore):
                if zmin_buff < component[k, i * n_ppore + j, 2] < zmax_buff:
                    dist = np.linalg.norm(component[k, i * n_ppore + j, :2] - pore_centers[:, i, k])
                    dist_from_center.append(dist)

    dist_from_center = np.array(dist_from_center)

    # Start setting up parameters necessary for binning
    bins = int(rmax/bin_width)  # the total number of bins
    r = np.linspace(0, rmax, int(bins) + 1)  # each discrete radius at which we will measure densities

    # Now bin the distances calculated above
    count = 0
    bin_contents_tails = np.zeros([int(bins) + 1])
    for i in range(len(dist_from_center)):
        distance = dist_from_center[i]
        if distance <= rmax:
            count += 1
            bin = int(np.floor((distance/rmax)*bins))
            bin_contents_tails[bin] += 1

    # calculate the density of ions in the area of the annulus defined by the inner and outer radii between which
    # components were calculated
    density = np.zeros([bins])
    for i in range(bins):
        density[i] = bin_contents_tails[i] / (math.pi*(((i + 1)*bin_width)**2 - (i*bin_width)**2))

    # normalize
    n = sum(density)
    for i in range(bins):
        density[i] /= n

    return density, r, bin_width


def gaus(x, a, sigma):
    return a*exp(-(x)**2/(2*sigma**2))

from scipy.misc import factorial


def poisson(k, lamb):
    return (lamb**k * np.exp(-lamb)) / factorial(k)

if __name__ == '__main__':

    args = initialize()  # parse the args

    t = md.load('%s' % args.input, top='%s' % args.gro)

    # Check for certain special arguments
    if args.component == 'tails':
        atoms = ['C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22',
                 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37',
                 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48']
    elif args.component == 'benzene':
        atoms = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
    elif args.component == 'tail_ends':
        atoms = ['C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C24', 'C25', 'C26', 'C27',
                 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45',
                 'C46', 'C47', 'C48']
    elif args.component == 'tail_fronts':
        atoms = ['C7', 'C8', 'C9', 'C21', 'C22', 'C23', 'C35', 'C36', 'C37']
    else:
        atoms = args.component

    if args.component != 'sys':
        atoms_to_keep = [a.index for a in t.topology.atoms if a.name in atoms]
        t.restrict_atoms(atoms_to_keep)
        pos = t.xyz
    else:
        # NOTE: if you use 'sys', use gmx trjconv -f *.trr -s *.tpr -pbc atom
        #                        then gmx trjconv -f *.trr -s *.tpr -pbc whole
        # This will make sure everything is in the box. The output is a .xtc file
        pos = t.xyz

    nT = np.shape(pos)[0]

    traj_start = int(args.start_frame)
    tot_atoms = np.shape(pos)[1]
    n_pores = int(args.pores)  # number of pores
    comp_ppore = tot_atoms/n_pores

    p_centers = avg_pore_loc(n_pores, pos, len(atoms))

    distances = 6  # number of p2p distances to calculate. My algorithm isn't smart enough for anything but six yet
    # distances = 16  # for 9 pore system
    p2ps = p2p(p_centers, distances)

    if args.auto_exclude:
        means = np.zeros([p2ps.shape[0]])
        for i in range(p2ps.shape[0]):
            means[i] = np.mean(p2ps[i, :])
        exclude = [np.argmax(means)]
    else:
        exclude = [int(i) for i in args.exclude]

    p2p_avg, p2p_std, equil = p2p_stats(p2ps, exclude, '%s' % args.nboot, '%s' % args.equil)
    print 'Average Pore to Pore distance: %.3f' % p2p_avg
    print 'Standard Deviation of Pore to Pore distances: %.3f' % p2p_std

    labels = ['1-2', '1-3', '1-4', '2-3', '2-4', '3-4']
    labels = ['1-4', '1-0', '0-2', '0-4', '0-3', '2-3', '2-5', '5-3', '3-4', '5-6', '3-7', '6-7', '4-7', '7-8',
              '4-8', '3-6']
    plt.figure(1)
    for i in range(distances):
        if i not in exclude:
            plt.plot(t.time, p2ps[i, :], label='%s' % labels[i])
    plt.title('Pore to Pore Distance Equilibration')
    plt.ylabel('Distance between pores (nm)')
    plt.xlabel('Time (ps)')
    # plt.legend(loc=1, fontsize=18)
    equil = 0
    density, r, bin_width = compdensity(pos, p_centers, equil, n_pores, buffer=0)

    for i in range(density.shape[-1]):
        if density[i] == 0:
            stop = i
            break

    # popt, pcov = curve_fit(gaus, r[:stop], density[:stop], p0=[density[0], 0.4])
    # # parameters, cov_matrix = curve_fit(poisson, r[:stop]/bin_width, density[:stop], p0=[1])
    # plt.figure(2)
    # plt.plot(r[:stop], gaus(r[:stop],*popt), 'ro:', label='fit')
    # # plt.plot(r[:stop], poisson(r[:stop]/bin_width, *parameters), 'r-', lw=2)
    # plt.title('Component Density Around Pore Center')
    # plt.xlabel('Distance from Pore Center (nm)')
    # plt.ylabel('Relative Component Density')
    # plt.bar(r[:stop], density[:stop], bin_width)
    if args.save:
        plt.savefig('p2p.png')
    if not args.noshow:
        plt.show()