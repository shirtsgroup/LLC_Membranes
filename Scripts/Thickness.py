#!/usr/bin/env python

# Find the membrane thickness based on a maximum and minimum z position
# Also gives the amount of space needed for water molecules given the amount of space wanted between periodic boundaries
# in the z direction

from __future__ import division
from __future__ import print_function
from past.utils import old_div
import argparse
from llclib import physical
import mdtraj as md
import matplotlib.pyplot as plt
from pymbar import timeseries
import numpy as np


def initialize():

    parser = argparse.ArgumentParser(description = 'Measure Membrane Thickness')

    parser.add_argument('-g', '--gro', type=str, default='wiggle.gro', help = 'Path to input file')
    parser.add_argument('-t', '--traj', type=str, default='traj_whole.xtc', help = 'Trajectory file (.xtc or .trr)')
    parser.add_argument('-w', '--water', default=6, help = 'nm of water between layers')
    parser.add_argument('-l', '--layers', default=20, help= 'Number of stacked layers in unit cell')
    parser.add_argument('-p', '--pores', default=4, help= 'Number of pores')
    parser.add_argument('-a', '--atoms', default=137, help='Number of atoms excluding ions')
    parser.add_argument('--trajectory', action="store_true")
    parser.add_argument('--ref_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms which will be used as'
                                                                 'a reference for determining max and min z coordinate')
    parser.add_argument('--save', action="store_true", help='Save plot output')
    parser.add_argument('--noshow', action="store_true", help='Do not display plot')
    parser.add_argument('-T', help='Use for concatenated trajectories of simulations run at multiple '
                                                 'temperatures. Supply the name of index file with temperatures and '
                                                 'times where the switch occurs. On each line supply the info with'
                                                 'the format temp:time')
    parser.add_argument('--plot_every', default=1, type=int, help='Plot every n frames')
    parser.add_argument('--plot_std', action="store_true", help='Plot average p2p distance at each frame')
    parser.add_argument('--nogrid', action="store_true", help='Do not calculate thickness based on a grid average')
    parser.add_argument('--gridres', default=1, type=int, help='Grid resolution. The box will be divided into'
                                                               'grid_res x grid_res equal sized boxes')
    parser.add_argument('--box', action="store_true", help='Calculate thickness based on box z dimension')

    args = parser.parse_args()

    return args


def thickness(filename):

    f = open(filename, "r")  # .gro file whose positions of Na ions will be read

    a = []  # list to hold lines of file
    for line in f:
        a.append(line)

    line = 0
    while a[line].count('HII') == 0:
        line += 1

    z = []  # list to hold z positions of all atoms

    benz_carbs = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']

    while a[line].count('HII') != 0:
        if str.strip(a[line][11:15]) in benz_carbs:
            z.append(float(a[line][36:44]))
        line += 1

    z_max = max(z)
    z_min = min(z)
    thick = z_max - z_min

    return thick, z_max, z_min


def parse_txt(txt, points, times, std=False):
    """
    Find out where to place lines in plot
    :param txt: a text file with pertinent information. Temperatures which the simulation was run at and at what time
    the simulation is at that temperature should be of the form Temperature:time (i.e. at `time` ps the simulation
    switches to `Temperature`). That section should be donoted by a header of the form [ T ]. You can also include a
    benchmark value which will be plotted as a straight line with dotted lines representing the error bounds. Start that
    section with the header [ bench ] with subsequent lines of the form value+/-error. Sections should be
    separated by a blank line
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
        y = np.array([0, 20])
        plt.plot(x, y, '--', linewidth=3, color='r', scaley=False)
        plt.annotate('T = %s' % T, xy=(time[i], np.mean(points)), xytext=(time[i] + 10, np.amax(points)))

    if std:
        std_equil = np.zeros([4, len(time) - 1])
        for i in range(len(time) - 1):
            begin = np.where(times == time[i])[0][0]
            end = np.where(times == time[i + 1])[0][0]
            slice = thick[begin:end]
            equil = timeseries.detectEquilibration(slice)[0]
            avg = np.mean(slice[equil:])
            std = np.std(slice[equil:])
            std_equil[:, i] = [avg, std, equil + begin, end]

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

    args = initialize()
    water_layer = int(args.water)  # nm of water wanted between layers

    if not args.trajectory:
        thick, z_max, z_min, thick_std = physical.thickness('%s' % args.gro, args.ref_atoms, not args.nogrid, grid_res=args.gridres)
        tot_thickness = thick + water_layer
        print('Membrane Thickness: %s +/- %s nm' % (thick, thick_std))
    else:

        # load trajectory and restrict atoms to reference atoms only
        t = md.load(args.traj, top=args.gro)
        keep = [a.index for a in t.topology.atoms if a.name in args.ref_atoms]
        pos = t.atom_slice(keep).xyz

        # Calculate a trajectory of thicknesses
        if args.box:
            thick = t.unitcell_vectors[:, 2, 2] #/ t.unitcell_vectors[:, 0, 0]
        else:
            thick = physical.thickness('%s' % args.gro, args.ref_atoms, not args.nogrid, pos)[0]

        # Decide when the system is equilibrated and take the average of thicknesses after equilibration
        equil_frame = timeseries.detectEquilibration(thick)[0]
        print('Equilibration detected after %d ns' % (old_div(t.time[equil_frame], 1000.0)))
        print('Average membrane thickness: %.2f +/- %.2f nm' % (np.mean(thick[equil_frame:]),
                                                                  np.std(thick[equil_frame:])))

        # Plot thickness vs time
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        # ax2 = ax1.twiny()
        ax1.plot(t.time[::args.plot_every]/1000, thick[::args.plot_every])
        #plt.title('Membrane thickness vs. time')
        ax1.set_xlabel('Time (ns)', fontsize=14)
        # new_tick_locations = np.linspace(0, t.time[-1], 7)
        # new_tick_labels = [int(280 + 60*(x/t.time[-1])) for x in new_tick_locations]
        # ax2.set_xlim(ax1.get_xlim())
        # ax2.set_xticks(new_tick_locations)
        # ax2.set_xticklabels(new_tick_labels)
        # ax2.set_xlabel("Temperature (K)", fontsize=14)
        ax1.set_ylabel('Membrane thickness (nm)', fontsize=14)
        ax1.tick_params(labelsize=14)
        # ax2.tick_params(labelsize=14)

        if args.T:
            xbounds, ybounds = parse_txt(args.T, thick[::args.plot_every], t.time[::args.plot_every]/1000, std=args.plot_std)
            plt.ylim(ybounds)
            plt.xlim(xbounds)

        # Save plot if you want
        if args.save:
            plt.savefig('thickness.png')

        # Don't show plot if you don't want to
        if not args.noshow:
            plt.show()
