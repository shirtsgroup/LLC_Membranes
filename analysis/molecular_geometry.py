#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
from scipy.spatial.distance import pdist


def initialize():

    parser = argparse.ArgumentParser(description='Model a continuous time random walk')

    # MD trajectory control
    parser.add_argument('-t', '--trajectory', default='PR_nojump.xtc', help='Path to input file.')
    parser.add_argument('-g', '--gro', default='em.gro', help='Name of .gro coordinate file.')

    return parser


class Molecule(object):

    def __init__(self, gro, traj):

        t = md.load(traj, top=gro)

        self.nframes = t.n_frames

        self.xyz = t.xyz

        self.radius = np.zeros([self.nframes])

    def radius(self):
        """ Calculate longest atom-atom distance at each frame of trajectory
        """

        for t in range(self.nframes):

            d = pdist(self.xyz[t, ...])
            print(d)
            exit()


if __name__ == "__main__":

    args = initialize()

    mol = Molecule(args.gro, args.trajectory)

    mol.radius()