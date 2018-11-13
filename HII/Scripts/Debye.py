#!/bin/bash

"""
    Apply the Debye scattering formula to atomic coordinates to produce a X-ray Diffraction pattern
"""

import argparse
from gro2xyz import coords_and_identity
import math
import Scattering_Factors
import numpy as np


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-f', '--file', default='NA_positions_5_images', help='Name of file to read')
    parser.add_argument('-o', '--output', default='wiggle.xyz', help='Name of file to be output')
    parser.add_argument('-l', '--lc', default='HII', help='Type of liquid crystal')
    parser.add_argument('-q', '--qrange', default=[0.05, 0.5], help='Range of q to be plotted')
    parser.add_argument('-s', '--step', default=0.01, help='Step size to take in q')

    args = parser.parse_args()
    return args


def scattering_factor(a, b, c, q):
    # d: 1/2d where d is the d-spacing
    sol = q/(4*math.pi)
    f = 0  # intialize summation
    for i in range(0, 4):  # for a1, b1, a2, b2, a3, b3, a4, b4
        f += a[i]*math.exp(-b[i]*sol**2)  # this function was used to fit the data in the Crystallographic tables
    f += c  # the constant c is contained outside the summation
    return f


def rij(coords, identity):

    atoms = np.shape(coords)[1]
    rij = np.zeros([atoms, atoms, atoms])

    for i in range(atoms):
        rij[i, , ]
        for k in range(atoms):



def diffraction_pattern(coords, identity, q):

    atoms = np.shape(coords)[1]

    I = np.zeros([len(q), 1])
    for k in range(len(q)):
        for i in range(atoms):
            exec "a = Scattering_Factors.%s.a" % identity[i]
            exec "b = Scattering_Factors.%s.b" % identity[i]
            exec "c = Scattering_Factors.%s.c" % identity[i]
            for j in range(atoms):
                f = scattering_factor(a, b, c, q[k])  # only one f since there is only one type of atom
                if rij[i, j] != 0:
                    I[k] += f*f*(math.sin(q[k]*rij[i, j])/q[k])
    return 0


if __name__ == "__main__":
    args = initialize()

    f = open('%s' % args.file, 'r')

    a = []
    for line in f:
        a.append(line)

    f.close()

    coords, identity = coords_and_identity(a, '%s' % args.lc)

    q_start = float(args.qrange[0])
    q_end = float(args.qrange[1])
    step = float(args.step)
    n_steps = (q_end - q_start) / step

    q = np.linspace(q_start, q_end, n_steps + 1)

    diffraction_pattern(coords, identity, q)