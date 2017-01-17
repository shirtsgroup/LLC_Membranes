#!/usr/bin/python

# This is the revamped version of Cylindricity.py

import numpy as np
import math
import matplotlib.pyplot as plt
from pylab import *
from scipy.interpolate import griddata
from matplotlib import animation
import argparse
import Get_Positions
from pymbar import timeseries
import random as ran


def initialize():
    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-i', '--input', default='wiggle_traj_bak.gro', help = 'Path to input file')
    parser.add_argument('-n', '--no_monomers', default=6, help = 'Number of Monomers per layer')
    parser.add_argument('-a', '--atoms', default=137, help = 'Number of atoms per monomer')
    parser.add_argument('-l', '--layers', default=20, help = 'Number of layers in each pore')
    parser.add_argument('-p', '--pores', default=4, help = 'Number of Pores')
    parser.add_argument('-c', '--component', default='sys', help = 'Counterion used to track pore positions')
    parser.add_argument('-f', '--start_frame', default=0, help = 'Frame number to start reading trajectory at')
    parser.add_argument('-s', '--layer_distribution', default='uniform', help = 'The distribution of monomers per layer')
    parser.add_argument('-L', '--alt_1', default=6, help = 'Monomers per layer for the first type of alternating layer')
    parser.add_argument('-A', '--alt_2', default=8, help = 'Monomers per layer for the second type of alternating layer')
    parser.add_argument('-d', '--direction', default='z', help = 'Axis along which to measure a component density')
    parser.add_argument('-S', '--slices', default='250', help = 'Number of slices to descretize chosen axis direction into')
    parser.add_argument('-g', '--grid_division', default=100, help = 'Number of blocks in x and y direction for heat map')
    parser.add_argument('-C', '--cmap', default='Blues', help = 'Color Scheme for heat map')
    parser.add_argument('-m', '--lc', default='HII', help = 'Type of liquid crystal monomer being used')
    parser.add_argument('-E', '--equil', default='auto', help = 'Frame number where system is equilibrated. "auto" will '
                        'use pymbar.timeseries.DetectEquilibration to determine which frame to start at. It is worth '
                        'double checking its choice manually')
    parser.add_argument('-x', '--exclude', default=4, help = 'Which pore-to-pore distance to exclude - pass the index of'
                                                            'the pore-to-pore distance as written in the list: '
                                                            '["1-2", "1-3", "1-4", "2-3", "2-4", "3-4"] ')
    parser.add_argument('-b', '--nboot', default=200, help = 'Number of bootstrap trials')

    args = parser.parse_args()
    return args


def avg_pore_loc(no_pores, nT, pos, comp_ppore):

    # Find the average location of the pores w.r.t. x and y

    p_center = np.zeros([2, no_pores, nT])

    for i in range(nT):
        for j in range(no_pores):
            for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                p_center[:, j, i] += pos[0:2, k, i]
            p_center[:, j, i] /= comp_ppore  # take the average

    return p_center


def p2p(p_centers, distances, nT):

        p2ps = np.zeros([distances, nT])  # distances in the order 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
        for i in range(nT):
            # So ugly ... sadness :(
            p2ps[0, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 1, i])
            p2ps[1, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 2, i])
            p2ps[2, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 3, i])
            p2ps[3, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 2, i])
            p2ps[4, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 3, i])
            p2ps[5, i] = np.linalg.norm(p_centers[:, 2, i] - p_centers[:, 3, i])

        return p2ps


def p2p_stats(p2ps, exclude, nboot, equil):

    frames = np.shape(p2ps)[1]
    exclude = int(exclude)
    nboot = int(nboot)

    # Get ride of the excluded trajectory
    p2p_new = np.zeros([5, frames])
    count = 0
    for i in range(6):
        if i != exclude:
            p2p_new[count, :] = p2ps[i, :]
            count += 1

    p2ps = p2p_new

    # Find the frame at which the system is equilibrated
    if equil == 'auto':
        ts = []
        for pore in range(5):
            ts.append(timeseries.detectEquilibration(p2ps[pore, :])[0])
        t = int(max(ts))  # use the max equil frame to ensure all pores are equilibrated
    else:
        t = int(equil)

    # Find the autocorrelation time for each pore - i.e. the time it takes for samples to become uncorrelated
    taus = []
    for i in range(5):
        tau = timeseries.integratedAutocorrelationTime(p2ps[i, t:])
        taus.append(tau)

    tau = int(max(taus))  # use the max again to ensure all trajectories are independent

    ind_trajectories = (frames - t) / tau  # the number of independent trajectories
    total_trajectories = ind_trajectories * 5  # the total number of trajectories
    trajectories = np.zeros([tau, total_trajectories])  # Create a new array to hold all the trajectories

    # fill up trajectory array
    for i in range(5):
        for j in range(ind_trajectories):
            trajectories[:, i * ind_trajectories + j] = p2ps[i, (t + j*tau):(t + (j + 1)*tau)]

    # bootstrap to get statistics
    p2p_boot = np.zeros([nboot])
    for i in range(nboot):
        p2p = 0
        for j in range(tau):
            T = ran.randrange(0, total_trajectories)  # pick a random trajectory from all the independent trajectories
            P = ran.randrange(0, tau)  # choose a random point in that trajectory
            p2p += trajectories[P, T]
        p2p_boot[i] = p2p / tau  # add the average from this trial to the p2p_boot

    # take average and standard deviation of bootstrap trial results
    p2p_avg = np.mean(p2p_boot)
    p2p_std = np.std(p2p_boot)

    return p2p_avg, p2p_std


if __name__ == '__main__':

    args = initialize()

    # p2ps = np.load('p2ps')
    # p2p_avg, p2p_std = p2p_stats(p2ps, '%s' % args.exclude, '%s' % args.nboot, '%s' % args.equil)
    # print p2p_avg, p2p_std
    # exit()
    # f = open(args.input, "r")  # .gro file whose positions of Na ions will be read
    # a = []  # list to hold lines of file
    # for line in f:
    #     a.append(line)
    # f.close()

    # pos, _, traj_pts, _ = Get_Positions.get_positions('%s' % args.input, '%s' % args.component, '%s' % args.lc,
    #                                                   'no')

    # atoms = ['C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22',
    #          'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37',
    #          'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48']
    #
    # # pos_tails, _, traj_pts, _ = Get_Positions.get_positions('%s' % args.input, atoms, '%s' % args.lc,
    # #                                                    'no')
    # # f = open('pos_tails', 'w')
    # # np.save(f, pos_tails)
    # # f.close()
    #
    # pos_tails = np.load('pos_tails')
    # pos = np.load('pos_NA')
    # pos_benzene = np.load('pos_benzene')
    # # traj_pts = np.shape(pos)[2]
    # # centers = np.zeros([3, 480, traj_pts])
    # # for k in range(traj_pts):
    # #     for i in range(480):
    # #         sum = np.array([0.0, 0.0, 0.0])
    # #         for j in range(6):
    # #             sum += pos[:, i*6 + j, k]
    # #         centers[:, i, k] = sum / 6.0
    #
    # # f = open('pos_benzene_centers', 'w')
    # # np.save(f, pos)
    # # f.close()
    # b_centers = np.load('pos_benzene_centers')
    # # p_centers = np.load('pcenters')
    # # print p_centers[:, :, 0]
    # # print b_centers[:, :6, 0]
    # # print np.linalg.norm(p_centers[:, 0, 0] - b_centers[:2, 0, 0])
    # # print np.linalg.norm(p_centers[:, 0, 0] - b_centers[:2, 118, 0])
    # # print np.linalg.norm(p_centers[:, 0, 0] - b_centers[:2, 117, 0])
    # # print np.linalg.norm(p_centers[:, 0, 0] - b_centers[:2, 116, 0])
    #
    #
    # # print np.shape(pos_tails)
    # # pos = np.load('pos_NA')
    # # pos_benzene = np.load('pos_benzene')
    # # pos_tails = np.load('pos_tails')
    # # print np.shape(pos_NA)
    # # print np.shape(pos_benzene)
    # # print np.shape(pos_tails)
    #
    # # pos = np.load('pos_array612ns')
    # # traj_pts = np.linspace(0, 614400, 1537)
    # traj_pts = np.shape(pos)[2]
    # # Questionably needed constants
    # no_atoms = int(args.atoms)  # number of atoms in a single monomer
    # no_layers = int(args.layers)  # number of layers in the membrane structure
    # mon_per_layer = int(args.no_monomers)  # number of monomers in each layer
    #
    # traj_start = int(args.start_frame)
    # tot_atoms = np.shape(pos)[1]
    # # nT = len(traj_pts)
    # nT = traj_pts
    # no_pores = int(args.pores)  # number of pores
    # # sim_length = traj_pts[-1]  # the last entry in the traj_pts list corresponds to the length of the simulation
    # comp_ppore = tot_atoms/no_pores
    #
    # # p_centers = avg_pore_loc(no_pores, nT, pos, comp_ppore)
    # # f = open('pcenters', 'w')
    # # np.save(f, p_centers)
    # p_centers = np.load('pcenters')
    # dist_from_center_NA = np.zeros([tot_atoms*traj_pts])
    # for k in range(nT):
    #     for i in range(4):
    #         for j in range(120):
    #             dist = np.linalg.norm(pos[:2, i * 120 + j, k] - p_centers[:, i, k])
    #             if dist < 2:
    #                 dist_from_center_NA[k * tot_atoms + i * 120 + j] = dist
    #             else:
    #                 dist_from_center_NA[k * tot_atoms + i * 120 + j] = 0.5
    #
    # tot_atoms = 2880
    # dist_from_center_benz = np.zeros([tot_atoms*traj_pts])
    # count = 0
    # for k in range(nT):
    #     for i in range(4):
    #         for j in range(720):
    #             dist = np.linalg.norm(pos_benzene[:2, i * 720 + j, k] - p_centers[:, i, k])
    #             dist_from_center_benz[k * tot_atoms + i * 720 + j] = dist
    #
    # # tot_atoms = 480 * len(atoms)
    # # n_ppore = 120 * len(atoms)
    # # print n_ppore
    # # dist_from_center_tails = np.zeros([tot_atoms*traj_pts])
    # # count = 0
    # # for k in range(nT):
    # #     for i in range(4):
    # #         for j in range(n_ppore):
    # #             dist = np.linalg.norm(pos_tails[:2, i * n_ppore + j, k] - p_centers[:, i, k])
    # #             dist_from_center_tails[k * tot_atoms + i * n_ppore + j] = dist
    # #
    # # f = open('all_tails', 'w')
    # # np.save(f, dist_from_center_tails)
    # # f.close()
    #
    # dist_from_center_tails = np.load('all_tails')
    # print np.shape(dist_from_center_tails)
    #
    max_x = 3.5
    bin_width = 0.1
    bins = int(max_x / bin_width)
    x_axis = np.linspace(bin_width/2, max_x - bin_width/2, int(bins))
    x_axis = np.linspace(0, max_x, int(bins) + 1)
    #
    # bin_contents_NA = np.zeros([bins])
    # for i in range(480 * 751):
    #     distance = dist_from_center_NA[i]
    #     if distance <= 3.5:
    #         bin = np.floor((distance / max_x)*bins)
    #         bin_contents_NA[bin] += 1
    #
    # NA_density = np.zeros([bins])
    # for i in range(bins):
    #     NA_density[i] = bin_contents_NA[i] / (math.pi*((i + 1)*bin_width)**2 - (i*bin_width)**2)
    #
    # bin_contents_benz = np.zeros([int(bins)])
    # for i in range(2880 * 751):
    #     distance = dist_from_center_benz[i]
    #     if distance <= 3.5:
    #         bin = int(np.floor((distance / max_x)*bins))
    #         bin_contents_benz[bin] += 1
    #
    # benz_density = np.zeros([bins])
    # for i in range(bins):
    #     benz_density[i] = bin_contents_benz[i] / (math.pi*((i + 1)*bin_width)**2 - (i*bin_width)**2)
    #
    # count = 0
    # bin_contents_tails = np.zeros([int(bins)])
    # for i in range(len(dist_from_center_tails)):
    #     distance = dist_from_center_tails[i]
    #     if distance <= 3.5:
    #         count += 1
    #         bin = int(np.floor((distance / max_x)*bins))
    #         bin_contents_tails[bin] += 1
    #
    # print count
    #
    # tails_density = np.zeros([bins])
    # for i in range(bins):
    #     tails_density[i] = bin_contents_tails[i] / (math.pi*((i + 1)*bin_width)**2 - (i*bin_width)**2)
    #
    # n_sodium = sum(NA_density)
    # n_benz = sum(benz_density)
    # n_tails = sum(tails_density)
    # for i in range(bins):
    #     tails_density[i] /= n_tails
    #     benz_density[i] /= n_benz
    #     NA_density[i] /= n_sodium
    upto = 25
    # f = open('tails_rho', 'w')
    # np.save(f, tails_density)
    # f.close()
    # f = open('benz_rho', 'w')
    # np.save(f, benz_density)
    # f.close()
    # f = open('NA_rho', 'w')
    # np.save(f, NA_density)
    # f.close()
    NA_density = np.load('NA_rho')
    benz_density = np.load('benz_rho')
    tails_density = np.load('tails_rho')

    plt.title('Component Density Around Pore Center', fontsize=26)
    plt.xlabel('Distance from Pore Center (nm)', fontsize=24)
    plt.ylabel('Relative Component Density', fontsize=24)
    font = {'family' : 'normal',
        'size'   : 16}
    plt.ylim(0, 0.28)
    plt.xlim(0, 2.5)
    matplotlib.rc('font', **font)
    plt.bar(x_axis[:upto], NA_density[:upto], bin_width, alpha=0.75, color='blue', label='Sodium Ions')
    plt.bar(x_axis[:upto], benz_density[:upto], bin_width, alpha=0.75, color='red', label='Benzene Carbons')
    plt.bar(x_axis[:upto], tails_density[:upto], bin_width, alpha=0.75, color=(0.25, 0.75, 0.75), label='Alkyl Chain Carbons')

    plt.legend(loc=1)
    plt.show()
    exit()

    weights_benz = np.ones_like(dist_from_center_benz) / len(dist_from_center_benz)
    # weights_NA = np.ones_like(NA_density) / len(NA_density)
    weights_tails = np.ones_like(dist_from_center_tails) / len(dist_from_center_tails)

    plt.title('Component Density around pore centers')
    # plt.hist(NA_density, bins = 50, alpha=0.33, label = 'NA', weights = weights_NA)
    plt.hist(dist_from_center_benz, bins = 50, alpha=0.33, label = 'Benzene Ring', weights=weights_benz)
    plt.hist(dist_from_center_tails, bins = 50, alpha=0.33, label = 'Alkyl Tails', weights=weights_tails)
    plt.xlabel('Distance from pore center (nm)')
    plt.ylabel('Fraction of group')
    plt.legend(loc=2)
    plt.show()
    exit()
    distances = 6  # number of p2p distances to calculate. My algorithm isn't good enough for anything but six yet
    p2ps = p2p(p_centers, distances, nT)
    # f = open('p2ps', 'w')
    # np.save(f, p2ps)
    # f.close()

    p2p_avg, p2p_std = p2p_stats(p2ps, '%s' % args.exclude, '%s' % args.nboot, '%s' % args.equil)

    labels = ['1-2', '1-3', '1-4', '2-3', '2-4', '3-4']
    plt.figure()
    for i in range(distances):
        plt.plot(traj_pts, p2ps[i, :], label='%s' % labels[i])

    plt.legend(loc=1, fontsize=18)
    plt.show()

