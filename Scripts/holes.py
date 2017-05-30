#!/usr/bin/env python

import numpy as np
from llclib import file_rw
import subprocess
import mdtraj as md
import argparse
import math
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Apply flat-bottomed position restraints based on a reference position')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='A coordinate file')

    args = parser.parse_args()

    return args


def position_restraints(file, atoms, axis, fconst):
    """
    Restrain the selected atoms in desired directions
    :param file: a list where each entry is a line from a coordinate file (.gro)
    :param atoms: indices of the atoms to position restrain
    :param axis: which direction to restrain
    :return: an array of position restraints formatted for easy writing into the topology (.itp)
    """

    # define force constants in their respective directions
    fcx = 0
    fcy = 0
    fcz = 0
    if 'x' in axis:
        fcx = fconst  # a large enough restraint to cause a large movement penalty
    if 'y' in axis:
        fcy = fconst
    if 'z' in axis:
        fcz = fconst

    atom_numbers = []  # find the numbers of the atoms which we are restraining
    with open(file, 'r') as f:
        for line in f:
            if str.strip(line[10:15]) in atoms:
                atom_numbers.append(int(line[15:20]))

    restraints = np.zeros([5, len(atom_numbers)])  # organize them into a list which can be translated to a topology
    for i in range(len(atom_numbers)):
        restraints[:, i] = [atom_numbers[i], 1, fcx, fcy, fcz]  # See: http://www.gromacs.org/Documentation/How-tos/Position_Restraints

    return restraints

if __name__ == "__main__":

    # box dimensions
    x = 8.0  # nm
    y = 8.0
    z = 8.0

    pore_radius = 0.5
    angle = np.pi / 3
    npts = 10000
    rows = 4
    npores = rows**2

    corners = np.zeros([8, 3])

    corners[0, :] = [0, 0, 0]
    corners[1, :] = [y*np.cos(angle), y*np.sin(angle), 0]
    corners[2, :] = [y*np.cos(angle) + x, y*np.sin(angle), 0]
    corners[3, :] = [y, 0, 0]
    for i in range(4, 8):
        corners[i, :] = corners[i - 4, :] + [0, 0, z]

    A = corners[0, :]
    AB = corners[1, :] - corners[0, :]  # v2 (referring to gromacs box formatting)
    AD = corners[3, :] - corners[0, :]  # v1
    AE = corners[4, :] - corners[0, :]  # v3

    pore_locations = np.zeros([npores, 2])
    p2p = x / rows

    for i in range(rows):
        for j in range(rows):
            pore_locations[i*rows + j, :] = [p2p*(0.5 + j) + p2p*(0.5 + i)*math.cos(angle), p2p/2*math.sin(angle) + i*p2p*math.sin(angle)]

    print pore_locations
    plt.scatter(pore_locations[:, 0], pore_locations[:, 1])
    plt.plot(corners[:4, 0], corners[:4, 1])
    plt.show()
    exit()
    points = np.zeros([npts, 3])

    # generate random points inside box
    for i in range(npts):
        u = np.random.rand()
        b = np.random.rand()
        h = np.random.rand()
        pt = A + u*AB + b*AD + h*AE
        points[i, :] = pt

    file_rw.write_gro_pos(points, 'holes.gro', box=[AD[0], AB[1], AE[2], AD[1], AD[2], AB[0], AB[2], AE[0], AE[1]])

    exit()

    pore_locations[0, :] = [x*.25, y*.25, 0]
    pore_locations[1, :] = [x*.75, y*.25, 0]
    pore_locations[2, :] = [x*.25, y*.75, 0]
    pore_locations[3, :] = [x*.75, y*.75, 0]

    n_ions = 10000  # ions to place in each pore
    positions = np.zeros([4*n_ions, 3])

    for i in range(4):
        for j in range(n_ions):
            zpos = np.random.rand()*z
            t = 2*np.pi*np.random.rand()
            u = np.random.rand()*pore_radius + np.random.rand()*pore_radius
            if u > pore_radius:
                r = 2*pore_radius - u
            else:
                r = u
            xpos = r*np.cos(t)
            ypos = r*np.sin(t)
            positions[i*n_ions + j, :] = [xpos, ypos, zpos] + pore_locations[i, :]

    file_rw.write_gro_pos(positions, 'holes.gro', box=[x, y, z])

    # dt = 0.0000001
    # p1 = subprocess.Poepn(["input.py", "-d", "%s" % dt, "-l", "%s" % (dt * 50), "--solvate", "-f", "50"])
    p = subprocess.Popen(["gmx", "solvate", "-cp", "holes.gro", "-cs", "spc216.gro", "-o", "solv.gro"])
    p.wait()

    t = md.load('solv.gro')
    keep = [a.index for a in t.topology.atoms if a.name != 'NA']
    water = t.atom_slice(keep)

    file_rw.write_gro(water, 'holes.gro')

    # apply restraints to the water

    # restraints = position_restraints('holes.gro', ['OW'], 'xy', 1000)
    #
    # with open('restraints.itp', 'w') as f:
    #     f.write("[ position_restraints ]\n")
    #     for i in range(restraints.shape[1]):
    #         f.write('{:6d}{:6d}{:1s}{:9f}{:1s}{:9f}{:1s}{:9f}\n'.format(int(restraints[0, i]), int(restraints[1, i]),'',
    #                                                         restraints[2, i], '', restraints[3, i], '', restraints[4, i]))