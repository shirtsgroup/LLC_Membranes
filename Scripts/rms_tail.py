#! /usr/bin/env python

import numpy as np
import argparse
import mdtraj as md
import lc_class


def initialize():

    parser = argparse.ArgumentParser(description='Write .mdp files and topology files for various types of simulation')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Name of trajectory')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-m', '--monomer', default='NAcarb11V.gro', type=str, help='Name of monomer coordinate file')
    parser.add_argument('-i', '--index', default='tail_index.ndx', help='plain text file containing names of all '
                                                                        'atoms of interest in each tail. See read_index'
                                                                        'function for format details')
    parser.add_argument('-a', '--axis', default='xy', help='Measure rmsd with respect to these axes')

    args = parser.parse_args()

    return args


def read_index(index):
    """
    :param index: index file
    Format: (1) Each group has a title formatted as follows: [ groupname ]
            (2) The line following each title should have the name of the atoms in each group written in order of their
            connectivity. They should all be contained on a single line
            (3) Each group should be separated by a blank line
    :return: index groups
    """

    grps = []
    sections = 0
    with open(index, 'r') as f:
        for line in f:
            if line.count('[ '):
                sections += 1
            elif line != "\n":
                grps.append(line.split())

    return grps


def lines(pos):
    """
    :param pos: xyz coordinates of all atoms in system (nframe, natoms, 3)
    :return: coefficients for equation of the line projecting out of the monitor for each monomer
    The equation is formatted parametrically of the form : [x, y, z] = [x0, y0, z0] + t[mx, my, mz]
    """

    nT = pos.shape[0]
    nlines = pos.shape[1] / 2
    coefficients = np.zeros([nT, nlines, 6])

    for t in range(nT):
        for i in range(nlines):
            m = pos[t, 2*i, :] - pos[t, 2*i + 1, :]
            coefficients[t, i, :3] = pos[t, 2*i + 1, :]
            coefficients[t, i, 3:6] = m

    return coefficients


def dist2line(x0, x1, x2, axis='xyz'):
    """
    See here for formula proof: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    :param x0: point anywhere in space. We are going to find out how far it is from the line
    :param x1: point on the line
    :param x2: another point on the line
    :return: minimum distance between x0 and line defined by x1 and x2
    """

    if axis == 'xy':
        x0 = x0[:2]
        x1 = x1[:2]
        x2 = x2[:2]

    a = x0 - x1
    b = x0 - x2
    c = x2 - x1
    axb = np.cross(a, b)

    return np.linalg.norm(axb) / np.linalg.norm(c)


def rmsd(deviations):
    """
    :param: deviations: an array of deviations
    :return: root mean square deviation of tail atoms from center line drawn through benzene
    """

    nT = deviations.shape[0]
    N = deviations.shape[1]
    rmsds = np.zeros([nT])

    for t in range(nT):
        sq = np.square(deviations[t, :])
        sum_sq = np.sum(sq)
        rmsds[t] = np.sqrt((1/float(N)) * sum_sq)

    return rmsds

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)  # load trajectory
    nT = t.n_frames

    mon = lc_class.LC(args.monomer)  # create monomer object
    lineatoms = mon.lineatoms  # indices of atoms in monomer used to make reference line. These happen to be the same
                               # ones that are used to build initial configurations
    nmon = t.xyz.shape[1] / mon.natoms  # number of monomers
    atoms_per_monomer = int(mon.natoms - mon.valence)  # atoms in each monomer excluding the counterion

    # isolate indices of line atoms in all monomers (2 per monomer)
    keep = [i for i in range(nmon*atoms_per_monomer) if i % atoms_per_monomer == lineatoms[0] or i % atoms_per_monomer == lineatoms[1]]
    pts = t.atom_slice(keep).xyz

    grps = read_index(args.index)
    atoms_per_tail = len(grps[0])

    tail1_index = [a.index for a in t.topology.atoms if a.name in grps[0]]
    tail3_index = [a.index for a in t.topology.atoms if a.name in grps[2]]

    # positions of all atoms in chains
    tail1 = t.atom_slice(tail1_index).xyz
    tail3 = t.atom_slice(tail3_index).xyz

    tail1_distance = np.zeros([nT, tail1.shape[1]])
    tail3_distance = np.zeros([nT, tail3.shape[1]])

    ntails = tail1.shape[1] / atoms_per_tail

    for t in range(nT):
        for i in range(ntails):
            for j in range(atoms_per_tail):
                tail1_distance[t, atoms_per_tail*i + j] = dist2line(tail1[t, atoms_per_tail*i + j], pts[t, 2*i, :], pts[t, 2*i + 1, :], axis=args.axis)
                tail3_distance[t, atoms_per_tail*i + j] = dist2line(tail3[t, atoms_per_tail*i + j], pts[t, 2*i, :], pts[t, 2*i + 1, :], axis=args.axis)

    RMSD = rmsd(tail1_distance)
    print RMSD