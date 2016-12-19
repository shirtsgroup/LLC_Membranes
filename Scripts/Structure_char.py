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
from pymbar.timeseries import detectEquilibration


def initialize():
    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-i', '--input', default='wiggle_traj_nopbc.gro', help = 'Path to input file')
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


if __name__ == '__main__':

    args = initialize()

    f = open(args.input, "r")  # .gro file whose positions of Na ions will be read
    a = []  # list to hold lines of file
    for line in f:
        a.append(line)
    f.close()

    # pos, _, traj_pts, _ = Get_Positions.get_positions('%s' % args.input, '%s' % args.component, '%s' % args.lc,
    #                                                   'no')
    pos = np.load('pos_array612ns')
    traj_pts = np.linspace(0, 614400, 1537)

    # Questionably needed constants
    no_atoms = int(args.atoms)  # number of atoms in a single monomer
    no_layers = int(args.layers)  # number of layers in the membrane structure
    mon_per_layer = int(args.no_monomers)  # number of monomers in each layer

    traj_start = int(args.start_frame)
    tot_atoms = np.shape(pos)[1]
    nT = len(traj_pts)
    no_pores = int(args.pores)  # number of pores
    sim_length = traj_pts[-1]  # the last entry in the traj_pts list corresponds to the length of the simulation
    comp_ppore = tot_atoms/no_pores

    p_centers = avg_pore_loc(no_pores, nT, pos, comp_ppore)

    print p_centers[:, :, -1]

    distances = 6  # number of p2p distances to calculate. My algorithm isn't good enough for anything but six yet
    p2ps = p2p(p_centers, distances, nT)
    f = open('p2ps', 'w')
    np.save(f, p2ps)
    f.close()

    print 'Calculating Equilibration Time'
    import time
    start = time.time()
    for i in range(0, 6):
        t, g, Neff_max = detectEquilibration(p2ps[i, :])
        print t
    print 'Equilibration detected in %s seconds' % (time.time() - start)

    labels = ['1-2', '1-3', '1-4', '2-3', '2-4', '3-4']
    plt.figure()
    for i in range(distances):
        plt.plot(traj_pts, p2ps[i, :], label='%s' % labels[i])

    plt.legend(loc=1, fontsize=18)
    plt.show()

