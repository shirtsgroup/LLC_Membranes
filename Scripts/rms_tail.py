#! /usr/bin/env python

import numpy as np
import argparse
import mdtraj as md
import lc_class
from pymbar import timeseries
import matplotlib.pyplot as plt
from tabulate import tabulate
import math
import tqdm
import progressbar


def initialize():

    parser = argparse.ArgumentParser(description='Write .mdp files and topology files for various types of simulation')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Name of trajectory')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-m', '--monomer', default='NAcarb11V.gro', type=str, help='Name of monomer coordinate file')
    parser.add_argument('-i', '--index', default='tail_index.ndx', help='plain text file containing names of all '
                                                                        'atoms of interest in each tail. See read_index'
                                                                        'function for format details')
    parser.add_argument('-a', '--axis', default='xy', help='Measure rmsd with respect to these axes')
    parser.add_argument('--time', action="store_true", help='analyze a time dependent trajectory')
    parser.add_argument('--taper', action="store_true", help='Step through each tail atom and calculate individual rmsds'
                                                             'then use that to calculate taper angle')
    parser.add_argument('--skip', default=1, type=int, help='Only read every nth frame')
    parser.add_argument('--avg', action="store_true", help='Calculate the average rmsd over all frames')
    parser.add_argument('--length', action="store_true", help='Calculate average length of monomers from head to tail')
    parser.add_argument('--all', action="store_true", help='Do all of the optional calculations')
    parser.add_argument('--noshow', action="store_true", help='Do not show plots')

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


def rotateplane(plane, pos, angle=0):
    """
    Calculate a rotation matrix to rotate a plane in 3 dimensions
    :param plane: indices of atoms making up plane which is being aligned in the xy plane
    :param pos: positions of atoms that need to be rotated
    :param angle: desired angle between xy plane (optional, default = 0 i.e. in plane)
    :return:
    """

    # vector pointing from point 1 to point 2
    v12 = plane[1, :] - plane[0, :]
    v13 = plane[2, :] - plane[0, :]

    # The cross product of v12 and v13 give a vector that is perpendicular to the plane:
    N = np.cross(v12, v13)
    N_desired = [0, math.sin(angle), math.cos(angle)]  # vector in the direction normal to our desired plane orientation

    RotationAxis = np.cross(N, N_desired)

    theta = math.acos(np.dot(N, N_desired)/(np.linalg.norm(N)*np.linalg.norm(N_desired)))  #  Rotation Angle (radians)

    L = [RotationAxis[0]/np.linalg.norm(RotationAxis), RotationAxis[1]/np.linalg.norm(RotationAxis),
                       RotationAxis[2]/np.linalg.norm(RotationAxis)]  # normalized Rotation Axis
    # ^ see: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/

    u, v, w = L[0], L[1], L[2]

    R = np.zeros((4, 4))
    R[3, 3] = 1
    R[0, 0] = u**2 + (v**2 + w**2)*math.cos(theta)  # math.cos takes theta in radians by default
    R[0, 1] = u*v*(1 - math.cos(theta)) - w*math.sin(theta)
    R[0, 2] = u*w*(1 - math.cos(theta)) + v*math.sin(theta)
    R[1, 0] = u*v*(1 - math.cos(theta)) + w*math.sin(theta)
    R[1, 1] = v**2 + (u**2 + w**2)*math.cos(theta)
    R[1, 2] = v*w*(1 - math.cos(theta)) - u*math.sin(theta)
    R[2, 0] = u*w*(1 - math.cos(theta)) - v*math.sin(theta)  # math.cos takes theta in radians by default
    R[2, 1] = w*v*(1 - math.cos(theta)) + u*math.sin(theta)
    R[2, 2] = w**2 + (u**2 + v**2)*math.cos(theta)

    b = np.ones([1])
    for i in range(pos.shape[0]):
        coord = np.concatenate((pos[i, :], b))
        x = np.dot(R, coord)
        pos[i, :] = x[:3]

    return pos


def dist2line(x0, x1, x2, plane, axis='xyz'):
    """
    See here for formula proof: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    Due to the way I calculate distance between the reference line and the tail carbon, the molecule needs to be rotated
    so that the benzene ring is coplanar with the xy axis. That way the xy projection is actually the xy distance b/w
    the reference line and atom.
    :param x0: point anywhere in space. We are going to find out how far it is from the line
    :param x1: point on the line
    :param x2: another point on the line
    :return: minimum distance between x0 and line defined by x1 and x2
    """

    if axis == 'xy':
        pts = np.zeros([3, 3])
        pts[0, :] = x0
        pts[1, :] = x1
        pts[2, :] = x2
        rot = rotateplane(plane, pts)
        x0r = rot[0, :2]
        x1r = rot[1, :2]
        x2r = rot[2, :2]
        a = x0r - x1r
        b = x0r - x2r
        c = x2r - x1r
    else:
        a = x0 - x1
        b = x0 - x2
        c = x2 - x1

    axb = np.cross(a, b)

    return np.linalg.norm(axb) / np.linalg.norm(c)


def pt2plane(x0, x1, x2, x3):
    """
    The plane is not rotated because the point to plane distance is exactly the distance above the plane of the monomer
    no matter how it is rotated, as long as the points defining the plane make up the benzene ring.
    :param x0: point whose distance from the plane we are interested in
    :param x1: known point in plane
    :param x2: known point in plane
    :param x3: known point in plane
    :return: distance of point x0 from plane defined by x1, x2 and x3
    """

    a = x2 - x1
    b = x3 - x1
    axb = np.cross(a, b)
    axb_norm = np.linalg.norm(axb)

    n = axb / axb_norm

    return np.dot(n, (x0 - x1))


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


def taper_line(tail, plane):
    """

    P1 and P2 are known points which define the reference line. We also know the position of all the carbons.
    To calculate the taper angle, find the angle defined by the two vectors representing the reference line and
    the path between C1 and C2 shown in the diagram below:
                                                 C2
                                                 ^
                                                 | d
      ref -- P1 --- P2 --------------------------v---

    :param tail: list of positions of all carbons on tail of interest [nT, natoms, 3]
    :return: slopes of line defined by first and last carbons of alkane chains with respect to reference line
    """

    tapers = np.zeros([nT, ntails])

    for t in range(nT):
        for i in range(ntails):
            # point on plane created by benzene. This plane will be rotated to be coplanar with xy plane
            plane = np.zeros([3, 3])
            plane[0, :] = planes[t, 3*i, :]
            plane[1, :] = planes[t, 3*i + 1, :]
            plane[2, :] = planes[t, 3*i + 2, :]
            # point which we need to rotate along with the plane
            points = np.zeros([4, 3])
            points[0, :] = pts[t, n*i, :]
            points[1, :] = pts[t, n*i + 1, :]
            points[2, :] = taper_ref[t, i, :]
            # points[2, :] = tail[t, atoms_per_tail * i, :]
            points[3, :] = tail[t, atoms_per_tail * (i + 1) - 1, :]
            # rotate points
            rot = rotateplane(plane, points)
            # extract projections onto xy plane
            pt1 = rot[0, :2]  # point in reference line
            pt2 = rot[1, :2]  # other point in reference line
            C1 = rot[2, :2]  # first carbon position
            C2 = rot[3, :2]  # last carbon position
            # calculate angle between vectors
            u = pt2 - pt1
            v = C2 - C1
            cos = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
            tapers[t, i] = np.arccos(cos) * (180 / np.pi)

    return tapers


def taper_plane(tail):
    """

    A reference plane is formed by the benzene ring. We are interested in the angle between a line formed between the
    first and last carbons on each chain, and the reference plane. This will give a 'z' taper angle.

    :param tail: list of positions of all carbons on tail of interest [nT, natoms, 3]
    :return: slopes of line defined by first and last carbons of alkane chains with respect to reference line
    """

    tapers = np.zeros([nT, ntails])

    for t in range(nT):
        for i in range(ntails):
            # find the normal vector to the plane
            # x1 = pts[t, n*i, :]
            # x2 = pts[t, n*i + 1, :]
            # x3 = pts[t, n*i + 2, :]
            # a = x2 - x1
            # b = x3 - x1
            # axb = np.cross(a, b)
            # axb_norm = np.linalg.norm(axb)
            # n = axb / axb_norm
            n = np.array([0, 0, 1])
            # find the vector in the direction pointing from the first carbon to the last carbon
            C1 = tail[t, atoms_per_tail * i, :]  # first carbon position
            C2 = tail[t, atoms_per_tail * (i + 1) - 1, :]  # last carbon position
            u = C2 - C1
            sin = abs(np.dot(n, u)) / (np.linalg.norm(n) * np.linalg.norm(u))
            tapers[t, i] = np.arcsin(sin) * (180 / np.pi)

    return tapers


def monomer_length(tail1, tail2, tail3, ref):
    """
    :param tail1: xyz positions of tail 1
    :param tail2: xyz positions of tail 2
    :param tail3: xyz positions of tail 3
    :param ref: reference atom position signifying the front end of the monomer
    :return: Average monomer length (nm)
    """

    nmon = ref.shape[1]  # number of monomers in unit cell
    atoms_per_tail = tail1.shape[1] / nmon  # number of atoms in each tail. Assuming identical length tails
    nT = ref.shape[0]  # number of frames

    d = []
    for t in range(nT):
        for i in range(nmon):
            ref_xyz = ref[t, i, :]
            # To speed up the calculation, only calculate distance between reference atom and last carbon on tail
            dist1 = np.linalg.norm(tail1[t, atoms_per_tail*(i + 1) - 1, :] - ref_xyz)
            dist2 = np.linalg.norm(tail2[t, atoms_per_tail*(i + 1) - 1, :] - ref_xyz)
            dist3 = np.linalg.norm(tail3[t, atoms_per_tail*(i + 1) - 1, :] - ref_xyz)
            d.append(dist1)
            d.append(dist2)
            d.append(dist3)

    return np.mean(d), np.std(d)

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[::args.skip]  # load trajectory
    time = t.time
    nT = t.n_frames
    bar = progressbar.ProgressBar()
    mon = lc_class.LC(args.monomer)  # create monomer object
    lineatoms = mon.lineatoms  # indices of atoms in monomer used to make reference line. These happen to be the same
                               # ones that are used to build initial configurations
    planeatoms = mon.planeatoms  # atoms which can be used to define the plane of the benzene ring
    taper_atom = mon.ref_atom_index

    nmon = t.xyz.shape[1] / mon.natoms  # number of monomers
    atoms_per_monomer = int(mon.natoms - mon.valence)  # atoms in each monomer excluding the counterion

    # isolate indices of line atoms in all monomers (2 per monomer)

    lines_keep = [i for i in range(nmon*atoms_per_monomer) if i % atoms_per_monomer == lineatoms[0] or i % atoms_per_monomer == lineatoms[1]]
    plane_keep = [a.index for a in t.topology.atoms if a.name in planeatoms]
    if args.taper or args.length or args.all:
        taper_keep = [i for i in range(nmon*atoms_per_monomer) if i % atoms_per_monomer == taper_atom]

    pts = t.atom_slice(lines_keep).xyz  # point defining reference line or plane
    planes = t.atom_slice(plane_keep).xyz  # point defining reference line or plane
    if args.taper or args.length or args.all:
        taper_ref = t.atom_slice(taper_keep).xyz  # points used as a reference for defining the taper angle

    grps = read_index(args.index)  # read the index file containing names of atoms in tail
    atoms_per_tail = len(grps[0])  # number of atoms in each tail (assumed to be equal in all tails)

    tail1_index = [a.index for a in t.topology.atoms if a.name in grps[0]]  # indices of all tail1 atoms
    if args.axis == 'z' or args.length or args.all:  # include the middle tail
        tail2_index = [a.index for a in t.topology.atoms if a.name in grps[1]]
    tail3_index = [a.index for a in t.topology.atoms if a.name in grps[2]]

    # positions of all atoms in chains
    tail1 = t.atom_slice(tail1_index).xyz
    if args.axis == 'z' or args.length or args.all:  # include the middle chain
        tail2 = t.atom_slice(tail2_index).xyz
    tail3 = t.atom_slice(tail3_index).xyz

    if args.axis == 'xy' or args.all:
        tail1_distance_xy = np.zeros([nT, tail1.shape[1]])
        tail3_distance_xy = np.zeros([nT, tail3.shape[1]])

    if args.axis == 'z' or args.all:
        tail1_distance_z = np.zeros([nT, tail1.shape[1]])
        tail2_distance_z = np.zeros([nT, tail2.shape[1]])
        tail3_distance_z = np.zeros([nT, tail3.shape[1]])

    ntails = tail1.shape[1] / atoms_per_tail  # really the number of monomers

    if args.axis == 'xy' or args.all:  # measure xy distance from each point in outside tails to reference line drawn down monomer
        print 'Calculating xy distance between carbon and reference line'
        for t in bar(range(nT)):
            for i in range(ntails):
                for j in range(atoms_per_tail):
                    plane = np.zeros([3, 3])
                    plane[0, :] = planes[t, 3*i, :]
                    plane[1, :] = planes[t, 3*i + 1, :]
                    plane[2, :] = planes[t, 3*i + 2, :]
                    tail1_distance_xy[t, atoms_per_tail*i + j] = dist2line(tail1[t, atoms_per_tail*i + j], pts[t, 2*i, :], pts[t, 2*i + 1, :], plane, axis=args.axis)
                    # tail2_distance[t, atoms_per_tail*i + j] = dist2line(tail2[t, atoms_per_tail*i + j], pts[t, 2*i, :], pts[t, 2*i + 1, :], axis=args.axis)
                    tail3_distance_xy[t, atoms_per_tail*i + j] = dist2line(tail3[t, atoms_per_tail*i + j], pts[t, 2*i, :], pts[t, 2*i + 1, :], plane, axis=args.axis)

            # mask arrays in case of nan's. They do not occur often but 1 is enough to mess up the whole calculation
            tail1_distance_xy = np.ma.masked_array(tail1_distance_xy, np.isnan(tail1_distance_xy))
            tail3_distance_xy = np.ma.masked_array(tail3_distance_xy, np.isnan(tail3_distance_xy))

    if args.axis == 'z' or args.all:  # measure distance from each tail to a plane defined by benzene ring
        print 'Calculating z distance between carbon and reference plane'
        for t in bar(range(nT)):
            for i in range(ntails):
                for j in range(atoms_per_tail):
                    tail1_distance_z[t, atoms_per_tail*i + j] = pt2plane(tail1[t, atoms_per_tail*i + j], planes[t, 3*i, :], planes[t, 3*i + 1, :], planes[t, 3*i + 2, :])
                    tail2_distance_z[t, atoms_per_tail*i + j] = pt2plane(tail2[t, atoms_per_tail*i + j], planes[t, 3*i, :], planes[t, 3*i + 1, :], planes[t, 3*i + 2, :])
                    tail3_distance_z[t, atoms_per_tail*i + j] = pt2plane(tail3[t, atoms_per_tail*i + j], planes[t, 3*i, :], planes[t, 3*i + 1, :], planes[t, 3*i + 2, :])

        # mask arrays in case of nan's. They do not occur often but 1 is enough to mess up the whole calculation
        tail1_distance_z = np.ma.masked_array(tail1_distance_z, np.isnan(tail1_distance_z))
        tail2_distance_z = np.ma.masked_array(tail2_distance_z, np.isnan(tail2_distance_z))
        tail3_distance_z = np.ma.masked_array(tail3_distance_z, np.isnan(tail3_distance_z))

    # if args.taper or args.all:
    #     RMSD1 = np.zeros([nT, atoms_per_tail])
    #     if args.axis == 'z' or args.all:
    #         RMSD2 = np.zeros([nT, atoms_per_tail])
    #     RMSD3 = np.zeros([nT, atoms_per_tail])
    #
    #     for i in range(atoms_per_tail):
    #         RMSD1[:, i] = rmsd(tail1_distance[:, i::atoms_per_tail])
    #         if args.axis == 'z' or args.all:
    #             RMSD2[:, i] = rmsd(tail2_distance[:, i::atoms_per_tail])
    #         RMSD3[:, i] = rmsd(tail3_distance[:, i::atoms_per_tail])
    #
    # else:
    #     RMSD1 = rmsd(tail1_distance)
    #     if args.axis == 'z':
    #         RMSD2 = rmsd(tail2_distance)
    #     RMSD3 = rmsd(tail3_distance)

    if args.length or args.all:

        L, L_std = monomer_length(tail1, tail2, tail3, taper_ref)

        print 'Average monomer length from head to tail : %.2f +/- %.2f' %(L, L_std)

    if args.taper or args.all:

        if args.axis == 'xy' or args.all:
            RMSD1_xy = np.zeros([nT, atoms_per_tail])
            RMSD3_xy = np.zeros([nT, atoms_per_tail])

        if args.axis == 'z' or args.all:
            RMSD1_z = np.zeros([nT, atoms_per_tail])
            RMSD2_z = np.zeros([nT, atoms_per_tail])
            RMSD3_z = np.zeros([nT, atoms_per_tail])

        for i in range(atoms_per_tail):

            if args.axis == 'xy' or args.all:
                RMSD1_xy[:, i] = rmsd(tail1_distance_xy[:, i::atoms_per_tail])
                RMSD3_xy[:, i] = rmsd(tail3_distance_xy[:, i::atoms_per_tail])

            if args.axis == 'z' or args.all:
                RMSD1_z[:, i] = rmsd(tail1_distance_z[:, i::atoms_per_tail])
                RMSD2_z[:, i] = rmsd(tail2_distance_z[:, i::atoms_per_tail])
                RMSD3_z[:, i] = rmsd(tail3_distance_z[:, i::atoms_per_tail])

        # find when each rmsd value equilibrates
        equil_frames = []

        for i in range(atoms_per_tail):  # look for equilibration of each tail atom

            if args.axis == 'xy' or args.all:
                equil_frames.append(timeseries.detectEquilibration(RMSD1_xy[:, i])[0])
                equil_frames.append(timeseries.detectEquilibration(RMSD3_xy[:, i])[0])

            if args.axis == 'z' or args.all:
                equil_frames.append(timeseries.detectEquilibration(RMSD1_z[:, i])[0])
                equil_frames.append(timeseries.detectEquilibration(RMSD2_z[:, i])[0])
                equil_frames.append(timeseries.detectEquilibration(RMSD3_z[:, i])[0])

        equil = max(equil_frames)  # choose the largest of them all to be safe

        if args.axis == 'xy' or args.all:
            RMSD_xy = np.zeros([atoms_per_tail])  # average all the values
            for i in range(atoms_per_tail):
                RMSD_xy[i] += np.mean(RMSD1_xy[equil:, i])
                RMSD_xy[i] += np.mean(RMSD3_xy[equil:, i])

            RMSD_xy /= 2.0

        if args.axis == 'z' or args.all:
            RMSD_z = np.zeros([atoms_per_tail])  # average all the values
            for i in range(atoms_per_tail):
                RMSD_z[i] += np.mean(RMSD1_z[equil:, i])
                RMSD_z[i] += np.mean(RMSD2_z[equil:, i])
                RMSD_z[i] += np.mean(RMSD3_z[equil:, i])

            RMSD_z /= 3.0

        # make a table
        headers = ['Carbon', 'RMSD']
        carbons = []
        for i in range(atoms_per_tail):
            carbons.append('C%s' % (i + 1))

        if args.axis == 'xy' or args.all:
            table = []
            if args.time:
                if nT > 20:  # not considering system equilibration since we are looking at a time dependent trajectory
                    jump = nT / 10
                    entries = 10
                else:
                    jump = 1
                    entries = nT
                for i in range(len(carbons)):
                    table.append([])
                    table[i].append(carbons[i])
                    for j in range(entries):
                        value = (np.mean(RMSD1_xy[j*jump:(j+1)*jump, i]) + np.mean(RMSD3_xy[j*jump:(j+1)*jump, i])) / 2
                        table[i].append(round(value, 3))
            else:
                for i in range(len(carbons)):
                    table.append([])
                    table[i].append(carbons[i])
                    table[i].append(round(RMSD_xy[i], 3))

            print 'RMSD of individual carbons with respect to xy axis'
            print tabulate(table, headers=headers)

        if args.axis == 'z' or args.all:
            table = []
            if args.time:
                if nT > 20:  # not considering system equilibration since we are looking at a time dependent trajectory
                    jump = nT / 10
                    entries = 10
                else:
                    jump = 1
                    entries = nT
                for i in range(len(carbons)):
                    table.append([])
                    table[i].append(carbons[i])
                    for j in range(entries):
                        value = (np.mean(RMSD1_z[j*jump:(j+1)*jump, i]) + np.mean(RMSD2_z[j*jump:(j+1)*jump, i]) +
                                 np.mean(RMSD3_z[j*jump:(j+1)*jump, i])) / 3
                        table[i].append(round(value, 3))
            else:
                for i in range(len(carbons)):
                    table.append([])
                    table[i].append(carbons[i])
                    table[i].append(round(RMSD_z[i], 3))

            print 'RMSD of individual carbons with respect to z axis'
            print tabulate(table, headers=headers)

        # now find the taper angle
        if args.axis == 'xy' or args.all:

            n = 2
            tapers_xy = np.zeros([nT, n*ntails])
            tapers_xy[:, :ntails] = taper_line(tail1, planes)
            tapers_xy[:, ntails:] = taper_line(tail3, planes)

        if args.axis == 'z' or args.all:

            n = 3
            tapers_z = np.zeros([nT, n*ntails])
            tapers_z[:, :ntails] = taper_plane(tail1)
            tapers_z[:, ntails:2*ntails] = taper_plane(tail2)
            tapers_z[:, 2*ntails:] = taper_plane(tail3)

        equil = []
        if args.axis == 'xy' or args.all:
            mean_xy = np.zeros([nT])
            for t in range(nT):
                mean_xy[t] = np.mean(tapers_xy[t, :])

            equil.append(timeseries.detectEquilibration(mean_xy)[0])

        if args.axis == 'z' or args.all:
            mean_z = np.zeros([nT])
            for t in range(nT):
                mean_z[t] = np.mean(tapers_z[t, :])

            equil.append(timeseries.detectEquilibration(mean_z)[0])

        equil = max(equil)

        if args.axis == 'xy' or args.all:
            print 'Average taper angle w.r.t in xy plane: %.2f degrees' % np.mean(tapers_xy[t, equil:])
            plt.plot(time, mean_xy)
        if args.axis == 'z' or args.all:
            print 'Average taper angle w.r.t z plane: %.2f degrees' % np.mean(tapers_z[t, equil:])
            plt.plot(time, mean_z)

        if not args.noshow:
            plt.show()

    if args.avg or args.all:

        if args.axis == 'xy' or args.all:

            RMSD1_xy = rmsd(tail1_distance_xy)
            RMSD3_xy = rmsd(tail3_distance_xy)

        if args.axis == 'z' or args.all:

            RMSD1_z = rmsd(tail1_distance_z)
            RMSD2_z = rmsd(tail2_distance_z)
            RMSD3_z = rmsd(tail3_distance_z)

        # need to do some modification of RMSD lists if you use --avg with --taper since the structure of the arrays are different : ([nT]) versus ([nT, mon_per_tail])
        equil_frames = []
        # print RMSD1.shape, RMSD2.shape, RMSD3.shape
        if args.axis == 'xy' or args.all:
            equil_frames.append(timeseries.detectEquilibration(RMSD1_xy)[0])
            equil_frames.append(timeseries.detectEquilibration(RMSD3_xy)[0])

        if args.axis == 'z' or args.all:
            equil_frames.append(timeseries.detectEquilibration(RMSD1_z)[0])
            equil_frames.append(timeseries.detectEquilibration(RMSD2_z)[0])
            equil_frames.append(timeseries.detectEquilibration(RMSD3_z)[0])

        equil = max(equil_frames)

        if args.axis == 'xy' or args.all:

            n = 2.0
            RMSDsum = RMSD1_xy[equil:] + RMSD3_xy[equil:]
            RMSD_xy = np.mean(RMSDsum) / n
            RMSD_xy_std = np.std(RMSDsum) / n
            print 'Root mean square deviation with respect to xy axis = %.3f +/- %.3f nm' % (RMSD_xy, RMSD_xy_std)

        if args.axis == 'z' or args.all:

            n = 3.0
            RMSDsum = RMSD1_z[equil:] + RMSD2_z[equil:] + RMSD3_z[equil:]
            RMSD_z = np.mean(RMSDsum) / n
            RMSD_z_std = np.std(RMSDsum) / n
            print 'Root mean square deviation with respect to z axis = %.3f +/- %.3f nm' % (RMSD_z, RMSD_z_std)

    if args.all:
        # special bonus unlockable when specifying the --all flag
        # The Volume of a monomer 'cone' with an elliptical base : (1/3)*pi*a*b*h
        # a : length  of semimajor axis (RMSD in xy)
        # b : length of semiminor axis (RMSD in z)
        # h : height of cone (monomer length)

        V = (1. / 3.) * np.pi * RMSD_z * RMSD_xy * L
        # uncertainty
        s = np.sqrt(L_std ** 2 + RMSD_z_std ** 2 + RMSD_xy_std ** 2)

        print 'Average monomer volume : %.3f +/- %.3f nm3' %(V, s)

    if args.time or args.all:

        if not args.avg or not args.all:

            if args.axis == 'xy':

                RMSD1_xy = rmsd(tail1_distance_xy)
                RMSD3_xy = rmsd(tail3_distance_xy)

            if args.axis == 'z':

                RMSD1_z = rmsd(tail1_distance_z)
                RMSD2_z = rmsd(tail2_distance_z)
                RMSD3_z = rmsd(tail3_distance_z)

        nfig = 0

        if args.axis == 'xy' or args.all:
            nfig += 1
            plt.figure(nfig)
            plt.plot(time, RMSD1_xy, label='Tail 1')
            plt.plot(time, RMSD3_xy, label='Tail 3')
            plt.title('Root mean square deviation of tail atoms from reference tail')
            plt.xlabel('time (ps)')
            plt.ylabel('RMSD')
            plt.legend()

        if args.axis == 'z' or args.all:
            nfig += 1
            plt.figure(nfig)
            plt.plot(time, RMSD1_z, label='Tail 1')
            plt.plot(time, RMSD2_z, label='Tail 2')
            plt.plot(time, RMSD3_z, label='Tail 3')
            plt.title('Root mean square deviation of tail atoms from reference plane')
            plt.xlabel('time (ps)')
            plt.ylabel('RMSD')
            plt.legend()

        if not args.noshow:
            plt.show()

