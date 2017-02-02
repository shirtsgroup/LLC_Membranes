#!/usr/bin/python

import numpy as np
import argparse


def initialize():

    parser = argparse.ArgumentParser(description='Find the distance between pores for a single frame')

    parser.add_argument('-i', '--input', default='wiggle.gro', type=str, help='Name of coordinate file to examine')

    args = parser.parse_args()

    return args

if __name__ == '__main__':

    args = initialize()

    f = open(args.input, 'r')

    a = []
    for line in f:
        a.append(line)

    f.close()

    pos = np.zeros([480, 3])

    count = 0
    for i in range(2, len(a) - 1):
        if str.strip(a[i][5:10]) == 'NA':
            pos[count, :] = [float(a[i][20:28]), float(a[i][28:36]), float(a[i][36:44])]
            count += 1

    centers = np.zeros([4, 2])

    for i in range(4):
        for j in range(120):
            centers[i, :] += pos[i*120 + j, :2]
        centers[i, :] /= 120
    print centers
    print np.linalg.norm(centers[0, :] - centers[1, :])
    print np.linalg.norm(centers[0, :] - centers[2, :])
    print np.linalg.norm(centers[0, :] - centers[3, :])
    print np.linalg.norm(centers[1, :] - centers[2, :])
    print np.linalg.norm(centers[1, :] - centers[3, :])
    print np.linalg.norm(centers[2, :] - centers[3, :])

