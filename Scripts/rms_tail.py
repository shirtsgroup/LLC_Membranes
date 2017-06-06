#! /usr/bin/env python

import numpy as np
import argparse
import mdtraj as md
import lc_class
from pymbar import timeseries
import matplotlib.pyplot as plt
from tabulate import tabulate
import math


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

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[::args.skip]  # load trajectory
    time = t.time
    nT = t.n_frames

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
    if args.taper:
        taper_keep = [i for i in range(nmon*atoms_per_monomer) if i % atoms_per_monomer == taper_atom]

    pts = t.atom_slice(lines_keep).xyz  # point defining reference line or plane
    planes = t.atom_slice(plane_keep).xyz  # point defining reference line or plane
    if args.taper:
        taper_ref = t.atom_slice(taper_keep).xyz  # points used as a reference for defining the taper angle

    grps = read_index(args.index)  # read the index file containing names of atoms in tail
    atoms_per_tail = len(grps[0])  # number of atoms in each tail (assumed to be equal in all tails)

    tail1_index = [a.index for a in t.topology.atoms if a.name in grps[0]]  # indices of all tail1 atoms
    if args.axis == 'z':  # include the middle tail
        tail2_index = [a.index for a in t.topology.atoms if a.name in grps[1]]
    tail3_index = [a.index for a in t.topology.atoms if a.name in grps[2]]

    # positions of all atoms in chains
    tail1 = t.atom_slice(tail1_index).xyz
    if args.axis == 'z':  # include the middle chain
        tail2 = t.atom_slice(tail2_index).xyz
    tail3 = t.atom_slice(tail3_index).xyz

    tail1_distance = np.zeros([nT, tail1.shape[1]])
    if args.axis == 'z': # include the middle chain
        tail2_distance = np.zeros([nT, tail2.shape[1]])
    tail3_distance = np.zeros([nT, tail3.shape[1]])

    ntails = tail1.shape[1] / atoms_per_tail  # really the number of monomers

    if args.axis == 'xy':  # measure xy distance from each point in outside tails to reference line drawn down monomer
        for t in range(nT):
            for i in range(ntails):
                for j in range(atoms_per_tail):
                    plane = np.zeros([3, 3])
                    plane[0, :] = planes[t, 3*i, :]
                    plane[1, :] = planes[t, 3*i + 1, :]
                    plane[2, :] = planes[t, 3*i + 2, :]
                    tail1_distance[t, atoms_per_tail*i + j] = dist2line(tail1[t, atoms_per_tail*i + j], pts[t, 2*i, :], pts[t, 2*i + 1, :], plane, axis=args.axis)
                    # tail2_distance[t, atoms_per_tail*i + j] = dist2line(tail2[t, atoms_per_tail*i + j], pts[t, 2*i, :], pts[t, 2*i + 1, :], axis=args.axis)
                    tail3_distance[t, atoms_per_tail*i + j] = dist2line(tail3[t, atoms_per_tail*i + j], pts[t, 2*i, :], pts[t, 2*i + 1, :], plane, axis=args.axis)

            # mask arrays in case of nan's. They do not occur often but 1 is enough to mess up the whole calculation
            tail1_distance = np.ma.masked_array(tail1_distance, np.isnan(tail1_distance))
            tail3_distance = np.ma.masked_array(tail3_distance, np.isnan(tail3_distance))

    elif args.axis == 'z':  # measure distance from each tail to a plane defined by benzene ring
        for t in range(nT):
            for i in range(ntails):
                for j in range(atoms_per_tail):
                    tail1_distance[t, atoms_per_tail*i + j] = pt2plane(tail1[t, atoms_per_tail*i + j], planes[t, 3*i, :], planes[t, 3*i + 1, :], planes[t, 3*i + 2, :])
                    tail2_distance[t, atoms_per_tail*i + j] = pt2plane(tail2[t, atoms_per_tail*i + j], planes[t, 3*i, :], planes[t, 3*i + 1, :], planes[t, 3*i + 2, :])
                    tail3_distance[t, atoms_per_tail*i + j] = pt2plane(tail3[t, atoms_per_tail*i + j], planes[t, 3*i, :], planes[t, 3*i + 1, :], planes[t, 3*i + 2, :])

        # mask arrays in case of nan's. They do not occur often but 1 is enough to mess up the whole calculation
        tail1_distance = np.ma.masked_array(tail1_distance, np.isnan(tail1_distance))
        tail2_distance = np.ma.masked_array(tail2_distance, np.isnan(tail2_distance))
        tail3_distance = np.ma.masked_array(tail3_distance, np.isnan(tail3_distance))

    if args.taper:
        RMSD1 = np.zeros([nT, atoms_per_tail])
        if args.axis == 'z':
            RMSD2 = np.zeros([nT, atoms_per_tail])
        RMSD3 = np.zeros([nT, atoms_per_tail])

        for i in range(atoms_per_tail):
            RMSD1[:, i] = rmsd(tail1_distance[:, i::atoms_per_tail])
            if args.axis == 'z':
                RMSD2[:, i] = rmsd(tail2_distance[:, i::atoms_per_tail])
            RMSD3[:, i] = rmsd(tail3_distance[:, i::atoms_per_tail])

    else:
        RMSD1 = rmsd(tail1_distance)
        if args.axis == 'z':
            RMSD2 = rmsd(tail2_distance)
        RMSD3 = rmsd(tail3_distance)

    if not args.time and not args.taper:

        equil_frames = []
        equil_frames.append(timeseries.detectEquilibration(RMSD1)[0])
        if args.axis == 'z':
            equil_frames.append(timeseries.detectEquilibration(RMSD2)[0])
        equil_frames.append(timeseries.detectEquilibration(RMSD3)[0])

        equil = max(equil_frames)

        RMSDsum = RMSD1[equil:] + RMSD3[equil:]
        if args.axis == 'z':
            RMSDsum += RMSD2[equil:]
            n = 3.0
        else:
            n = 2.0

        RMSD = np.mean(RMSDsum) / n

        print 'Root mean square deviation with respect to %s axis = %.3f nm' % (args.axis, RMSD)

    elif args.taper:

        # find when each rmsd value equilibrates
        equil_frames = []

        for i in range(atoms_per_tail):  # look for equilibration of each tail atom
            equil_frames.append(timeseries.detectEquilibration(RMSD1[:, i])[0])
            if args.axis == 'z':
                equil_frames.append(timeseries.detectEquilibration(RMSD2[:, i])[0])
            equil_frames.append(timeseries.detectEquilibration(RMSD3[:, i])[0])

        equil = max(equil_frames)  # choose the largest of them all to be safe

        RMSD = np.zeros([atoms_per_tail])  # average all the values
        for i in range(atoms_per_tail):
            RMSD[i] += np.mean(RMSD1[equil:, i])
            if args.axis == 'z':
                RMSD[i] += np.mean(RMSD2[equil:, i])
            RMSD[i] += np.mean(RMSD3[equil:, i])

        if args.axis == 'z':
            n = 3  # 3 tails are analyzed when looking at the z axis
        else:
            n = 2

        RMSD /= n

        # make a table
        headers = ['Carbon', 'RMSD']
        carbons = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10', 'C11']

        table = []
        for i in range(len(carbons)):
            table.append([])
            table[i].append(carbons[i])
            table[i].append(round(RMSD[i], 3))

        print tabulate(table, headers=headers)

        # now find the taper angle
        if args.axis == 'xy':

            tapers = np.zeros([nT, 2*ntails])
            tapers[:, :ntails] = taper_line(tail1, planes)
            tapers[:, ntails:] = taper_line(tail3, planes)

        if args.axis == 'z':

            tapers = np.zeros([nT, 3*ntails])
            tapers[:, :ntails] = taper_plane(tail1)
            tapers[:, ntails:2*ntails] = taper_plane(tail2)
            tapers[:, 2*ntails:] = taper_plane(tail3)

        mean = np.zeros([nT])
        for t in range(nT):
            mean[t] = np.mean(tapers[t, :])

        equil = timeseries.detectEquilibration(mean)[0]

        print 'Average taper angle: %.2f degrees' % np.mean(tapers[t, equil:])
        plt.plot(time, mean)
        plt.show()

    if args.time:

        if args.taper:
            RMSD1 = rmsd(tail1_distance)
            if args.axis == 'z':
                RMSD2 = rmsd(tail2_distance)
            RMSD3 = rmsd(tail3_distance)

        plt.plot(time, RMSD1, label='Tail 1')
        if args.axis == 'z':
            plt.plot(time, RMSD2, label='Tail 2')
            reference = 'plane'
        else:
            reference = 'tail'
        plt.plot(time, RMSD3, label='Tail 3')
        plt.title('Root mean square deviation of tail atoms from reference %s' % reference)
        plt.xlabel('time (ps)')
        plt.ylabel('RMSD')
        plt.legend()
        plt.show()

