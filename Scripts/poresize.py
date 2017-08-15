#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
from past.utils import old_div
import argparse
import mdtraj as md
import numpy as np
from llclib import physical
import matplotlib.pyplot as plt
from pymbar import timeseries


def initialize():

    parser = argparse.ArgumentParser(description='Measure the density inside the hydrophobic region of an LLC system')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help = 'Trajectory file (.xtc or .trr). The'
                        'trajectory should be preprocessed beforehand to make molecules whole or else the distance '
                        'calculation get thrown off. Use "gmx trjconv -pbc whole" to do this')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help = 'Name of coordinate file')
    parser.add_argument('-c', '--components', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'],
                        help='component used to track pore positions and define pore radius')
    parser.add_argument('--noshow', action="store_true", help='Do not show the plots')
    parser.add_argument('--save', action="store_true", help='Save plots')
    parser.add_argument('--normalize', action="store_true", help='Normalize order parameter so max is 1 and min is 0')
    parser.add_argument('-T', help='Use for concatenated trajectories of simulations run at multiple '
                                                 'temperatures. Supply the name of index file with temperatures and '
                                                 'times where the switch occurs. On each line supply the info with'
                                                 'the format temp:time')
    parser.add_argument('--plot_every', default=1, type=int, help='Plot every n frames')
    parser.add_argument('--plot_std', action="store_true", help='Plot average p2p distance at each frame')

    args = parser.parse_args()

    return args


def parse_txt(txt, points, times, name, std=False):
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
        plt.plot(x, y, '--', linewidth=3, color='r')
        plt.annotate('T = %s' % T, xy=(time[i], np.mean(points)), xytext=(time[i] + 10000, np.amax(points)))

    if std:
        std_equil = np.zeros([4, len(time) - 1])
        for i in range(len(time) - 1):
            begin = np.where(times == time[i])[0][0]
            end = np.where(times == time[i + 1])[0][0]
            slice = points[begin:end]
            equil = timeseries.detectEquilibration(slice)[0]
            avg = np.mean(slice[equil:])
            std = np.std(slice[equil:])
            print(avg, std)
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
    while a[B].count('[ bench_%s ]' % name) == 0:
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
        plt.plot(x, y, linewidth=2, color='g')
        # plt.plot(x, upperbound, '--', linewidth=1, color='b')
        # plt.plot(x, lowerbound, '--', linewidth=1, color='b')

    ybounds = np.array([min(y_low)*.99, max(y_high)*1.01])
    xbounds = np.array([-100, max(times)])

    return xbounds, ybounds

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)  # load gromacs trajectory
    keep = [a.index for a in t.topology.atoms if a.name in args.components]   # get the index of all atoms in components
    pore_components = t.atom_slice(keep)  # create a new trajectory object only describing those atoms
    pos = pore_components.xyz

    pcenters = physical.avg_pore_loc(4, pos)  # find the pore centers

    r, r_std = physical.limits(pos, pcenters)

    poresize_equil = timeseries.detectEquilibration(r)[0]
    order_equil = timeseries.detectEquilibration(old_div(r,r_std))[0]

    if args.normalize:
        order = old_div(r, r_std)
        order -= min(order)
        order /= max(order)

    print('Pore size equilibrated after %d ns' % (old_div(t.time[poresize_equil], 1000)))
    print('Average Pore Size: %.2f +/- %.2f nm' %(np.mean(r[poresize_equil:]), np.mean(r_std[poresize_equil:])))
    print('Order parameter equilibrated after %d ns' % (old_div(t.time[order_equil], 1000)))
    if args.normalize:
        print('Average Order Parameter: %.2f' % np.mean(order[order_equil:]))
    else:
        print('Average Order Parameter: %.2f +/- %.2f' % (np.mean(old_div(r[order_equil:],r_std[order_equil:])),
                                                          np.std(old_div(r[order_equil:],r_std[order_equil:]))))

    plt.figure(1)
    plt.plot(t.time[::args.plot_every], r[::args.plot_every])
    plt.title('Pore size vs time')
    plt.ylabel('Average distance of benzene from pore center (nm)')
    plt.xlabel('Time')

    if args.T:
        xbounds, ybounds = parse_txt(args.T, r[::args.plot_every], t.time[::args.plot_every], 'poresize', std=args.plot_std)
        plt.ylim(ybounds)
        plt.xlim(xbounds)

    if args.save:
        plt.savefig('poresize.png')

    # plt.figure(2)
    # plt.plot(t.time, r_std)
    # plt.ylabel('Standard deviation of distance of benzene from pore center')
    # plt.xlabel('Time')

    plt.figure(3)

    T = np.linspace(280, 340, t.time.shape[0])

    if args.normalize:
        # plt.plot(t.time, order)
        plt.plot(T, order)
    else:
        plt.plot(t.time[::args.plot_every], old_div(r,r_std)[::args.plot_every])

    plt.title('Order parameter vs time')
    plt.xlabel('Time (ps)')
    # plt.xlabel('Temperature (K)')
    plt.ylabel('Order parameter')

    if args.T:
        xbounds, ybounds = parse_txt(args.T, old_div(r,r_std)[::args.plot_every], t.time[::args.plot_every], 'order', std=args.plot_std)
        plt.ylim(ybounds)
        plt.xlim(xbounds)

    if args.save:
        plt.savefig('order.png')

    if not args.noshow:
        plt.show()