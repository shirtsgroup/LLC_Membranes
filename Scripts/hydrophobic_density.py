#! /usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
from llclib import physical, file_rw
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Measure the density inside the hydrophobic region of an LLC system')

    parser.add_argument('-t', '--traj', default='wiggle.xtc', type=str, help = 'Trajectory file (.xtc or .trr). The'
                        'trajectory should be preprocessed beforehand to make molecules whole or else the distance '
                        'calculation get thrown off. Use "gmx trjconv -pbc whole" to do this')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help = 'Name of coordinate file')
    parser.add_argument('-i', '--inner_comps', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'],
                        help = 'component used to track pore positions and define inner limit of alkane region')
    parser.add_argument('-o', '--outer_comps', nargs='+', default=['C17', 'C31', 'C45'], help = 'Name of tail components '
                                                                        'used to define outer limit of alkane region')
    parser.add_argument('--vis', help='visualize the atoms that are included in the density calculation by outputting'
                                      'a restricted .gro file', action="store_true")
    parser.add_argument('--noshow', help='specify this flag if you do not want a plot shown after the calculations. The'
                                         'average of the density will then be reported', action="store_true")
    parser.add_argument('--nudge_inner', default = 1, type=float, help='number of sigma to move the inner limit out')
    parser.add_argument('--nudge_outer', default= 1, type=float, help='number of sigma to move the outer limit in')

    args = parser.parse_args()

    return args


def limits(pos, pcenters):
    """
    Estimate the pore 'radius' based on the position of some component and it's maximum deviation from the pore center
    :param: pos: the positions of all atoms included in making the estimate
    :param: pcenters: the x,y positions of the pore centers for each frame
    :return: an approximate pore radius. Beyond which, we have entered the alkane region
    """

    nT = pcenters.shape[2]
    npores = pcenters.shape[1]
    natoms= pos.shape[1]
    atom_ppore = natoms / npores

    deviation = np.zeros([nT, npores, atom_ppore])
    for f in range(nT):
        for i in range(atom_ppore):
            for j in range(npores):
                deviation[f, j, i] = np.linalg.norm(pos[f, j*atom_ppore + i, :2] - pcenters[:, j, f])  # left off here

    deviation = np.reshape(deviation, (nT, natoms))
    fr = np.zeros([nT])
    frstd = np.zeros([nT])

    for i in range(nT):
        fr[i] = np.mean(deviation[i, :])  # + np.std(deviation[i, :]) # maybe?
        frstd = np.std(deviation[i, :])

    return fr, frstd


def density(inner_limits, outer_limits, pcenters, t):
    """
    :param inner_limits: The inner limits to the alkane region for each frame
    :param outer_limits: The outer limits to the alkane region for each frame
    :param t: An mdtraj trajectory object for the system
    :param pcenters: the location of the pore centers at each frame
    :return: density of alkane region at each frame
    """

    NA = 6.022*10**23  # avogadros number
    conv = 1*10**-21  # convert nm^3 to cm^3
    pos = t.xyz  # get positions
    box = t.unitcell_vectors


    nT = pos.shape[0]
    atoms = pos.shape[1]
    npores = pcenters.shape[1]
    atomsppore = atoms / npores

    # calculate volume of alkane region at each frame
    vol = np.zeros([nT])
    for i in range(nT):
        a = np.pi*(outer_limits[i]**2 - inner_limits[i]**2) * npores # cross sectional area of annulus defined by limits
        # if there is worry about overlap of the circles, here is an algorithm that can be implemented:
        # http://stackoverflow.com/questions/1667310/combined-area-of-overlapping-circles
        z = np.linalg.norm(box[i, 2, :])  # z dimension

        vol[i] = a*z*conv

    mass = np.zeros([atoms])  # molar mass of each atom
    count = 0
    for a in t.topology.atoms:
        mass[count] = a.element.mass
        count += 1

    d = np.zeros([nT])
    for f in range(nT):
        for j in range(npores):
            for k in range(atomsppore):
                if inner_limits[f] < np.linalg.norm(pos[f, j*atomsppore + k, :2] - pcenters[:, j, f]) \
                        < outer_limits[f]:  # check if x^2 + y^2 is inside limits
                    d[f] += (mass[j*atomsppore + k] / NA)  # add the mass of the atom meeting conditions and convert to mass of single atom

    for f in range(nT):
        d[f] /= vol[f]

    return d


def visualize(inner_limits, outer_limits, pcenters, t, frame=-1):
    """
    Create a .gro file so we can see which atoms are included in the density calculation (single frame)
    :param inner_limits: The inner limits to the alkane region for each frame
    :param outer_limits: The outer limits to the alkane region for each frame
    :param pcenters: the location of the pore centers at each frame
    :param t: An mdtraj trajectory object for the system
    :param frame: frame number you want to look at, default is last frame
    :return: a new trajectory object only containing atoms that are within the defined limits
    """

    pos = t.xyz
    natoms = pos.shape[1]
    npores = pcenters.shape[1]
    atomsppore = natoms / npores

    keep = []
    for i in range(npores):
        for j in range(atomsppore):
            if inner_limits[frame] < np.linalg.norm(pos[frame, i*atomsppore + j, :2] - pcenters[:, i, frame]) \
                    < outer_limits[frame]:
                keep.append(i*atomsppore + j)

    return t.atom_slice(keep)

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)  # load trajectory using mdtraj

    # get indices of component that will be use to find pore centers and define inner limit of alkane region
    innerlist = [a.index for a in t.topology.atoms if a.name in args.inner_comps]
    outerlist = [a.index for a in t.topology.atoms if a.name in args.outer_comps]

    inner = t.atom_slice(innerlist).xyz
    outer = t.atom_slice(outerlist).xyz

    p_centers = physical.avg_pore_loc(4, inner)  # these will become the centers of the circle I create to define the pore region

    inner_limits, inner_std = limits(inner, p_centers)
    outer_limits, outer_std = limits(outer, p_centers)

    inner_limits += inner_std * args.nudge_inner
    outer_limits -= outer_std * args.nudge_outer

    d = density(inner_limits, outer_limits, p_centers, t)

    if args.vis:
        restricted = visualize(inner_limits, outer_limits, p_centers, t)
        file_rw.write_gro(restricted, 'restricted.gro')

    if not args.noshow:
        plt.plot(t.time, d)
        plt.show()
    else:
        print "Average Density: %s g/cm^3" % np.mean(d[-5:])

