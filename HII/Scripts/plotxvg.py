#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
from pymbar.timeseries import detectEquilibration


def initialize():

    parser = argparse.ArgumentParser(description='Plot data in .xvg file')

    parser.add_argument('-f', '--xvg', default='energy.xvg', type=str, help='Output file from gmx energy')
    parser.add_argument('-u', '--units', default='ps', type=str, help='Units to plot time axis')
    parser.add_argument('-s', '--save', action="store_true", help='Save plot')
    parser.add_argument('-o', '--out', default='xvg.png', type=str, help='Name to save plot under')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    xvg = []
    with open(args.xvg, 'r') as f:
        for line in f:
            xvg.append(line)

    wham = False
    if xvg[2].count('wham') != 0:
        wham = True

    if wham:
        start_string = '@TYPE'
    else:
        start_string = '@ s'

    start = 0
    while xvg[start].count('@') == 0:
        start += 1

    title = xvg[start].split('"')[1]
    xlabel = xvg[start + 1].split('"')[1]
    ylabel = xvg[start + 2].split('"')[1]

    while xvg[start].count(start_string) == 0:
        start += 1

    if not wham:
        columns_labels = [xlabel]
        while xvg[start].count(start_string) == 1:
            columns_labels.append(xvg[start].split('"')[1])
        start += 1

    if wham:
        start += 1
        columns_labels = []
        for i in range(len(xvg[start].split())):
            columns_labels.append('Column %s' % (i + 1))

    data = np.zeros([len(xvg) - start, len(columns_labels)])
    for i in range(start, len(xvg)):
        linedata = [float(d) for d in xvg[i].split()]
        data[i - start, :] = linedata

    equil = np.zeros([data.shape[1] - 1])  # equilibration time for each timeseries
    for i in range(1, data.shape[1]):
        equil[i - 1] = detectEquilibration(data[:, i])[0]

    if args.units == 'ns':
        data[:, 0] /= 1000

    plt.figure()
    for i in range(1, len(columns_labels)):
        plt.plot(data[:, 0], data[:, i], label=columns_labels[i], linewidth=2)
    if not wham:
        plt.legend()
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    if args.save:
        plt.savefig(args.out)

    for i in range(equil.size):
        print('Average %s: %.2f +\- %.2f' % (columns_labels[i + 1], np.mean(data[int(equil[i]):, i + 1]), np.std(data[int(equil[i]):, i + 1])))

    plt.show()