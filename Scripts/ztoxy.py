#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import subprocess
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Grab the last trajectory frame and make a .gro file')

    parser.add_argument('-e', '--edr', default='wiggle.edr', type=str, help='Name of gromacs energy file')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    box_pipes = [20, 21]
    x = []
    y = []

    for i in range(len(box_pipes)):

        x.append([])
        y.append([])

        ps = subprocess.Popen(('echo', '%s' % box_pipes[i]), stdout=subprocess.PIPE)
        subprocess.call(['gmx', 'energy', '-f', '%s' % args.edr], stdin=ps.stdout)
        ps.wait()

        a = []
        with open('energy.xvg', 'r') as f:

            for line in f:
                a.append(line)

        begindex = 0
        while a[begindex].count('#') == 1 or a[begindex].count('@') == 1:
            begindex += 1

        for j in range(begindex, len(a)):
            data = a[j].split()
            x[i].append(data[0])
            y[i].append(data[1])

    ratio = [float(y[0][i]) / float(y[1][i]) for i in range(len(y[0]))]

    plt.plot(x[0], ratio)
    plt.show()