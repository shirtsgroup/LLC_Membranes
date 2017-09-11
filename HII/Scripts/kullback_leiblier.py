#!/usr/bin/env python

import argparse
import math
import numpy as np
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Edit topology of a system or residue to repartition mass')

    parser.add_argument('-p', '--dist1', type=str, help='Distribution stored in a pickled numpy array')
    parser.add_argument('-q', '--dist2', type=str, default='uniform', help='Second distribution for comparison. The default, '
                                                                 'uniform, will compare to normal distribution')
    args = parser.parse_args()

    return args


def kullback_leiblier(p, q):
    """
    Find the kullback-leiblier divergence from q to p

    :param p: probability distribution 1
    :param q: probability distribution 2
    :return: kb : kullback-leiblier divergence - the amount of information lost when q is used to approximate p
    """
    # normalize
    p /= sum(p)
    q /= sum(q)

    nbins = p.shape[0]

    kb = 0
    for i in range(nbins):
        if p[i] != 0 and q[i] != 0:
            kb += p[i]*math.log(p[i]/q[i])

    return kb

if __name__ == "__main__":

    args = initialize()

    dist = np.load(args.dist1)
    angles = dist['angles']
    nbins = 45
    angles = [value for value in angles if not math.isnan(value)]
    (p, bins, patches) = plt.hist(angles, bins=nbins)

    if args.dist2 == 'uniform':

        q = np.zeros_like(p) + 1

    print(kullback_leiblier(p, q))
