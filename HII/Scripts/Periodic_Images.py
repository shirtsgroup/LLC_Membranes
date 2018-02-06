#!/usr/bin/env python

""" Use this script to duplicate points periodically in the x-y directions.

    Example: images = 1, duplicate and shift 'point' into each number grid
  ___________________________
  \        \        \        \
   \   I    \   II   \  III   \
    \________\________\________\
     \        \        \        \
      \  IV    \  point \   V    \
       \________\________\________\
        \        \        \        \
         \   VI   \  VII   \  VIII  \
          \________\________\________\

    images = 1 corresponds to the above 3 x 3 grid. images = 2 corresponds to a 5 x 5 grid. images = 3 is 7 x 7 etc.
    The meaning of images is the number of times we completely surround point with duplicates of itself

    x_box and y_box are defined as follows:
    <-x_box->
    __________  ^
    \        \   \
     \  box   \  y_box
      \________\  \
                   v
"""
from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
import argparse
import numpy as np
import mdtraj as md
from llclib import file_rw


def initialize():

    parser = argparse.ArgumentParser(description='Duplicate points periodically in the x-y directions')

    parser.add_argument('-f', '--file', default='15nm.gro', help='File to replicate periodically')
    parser.add_argument('-o', '--output', default='Periodic.gro', help='Name of output file')
    parser.add_argument('-i', '--images', default=1, help='Number of periodic images')
    parser.add_argument('-a', '--angle', default=60, help='Angle between x and y box vector')

    args = parser.parse_args()

    return args


def shift_matrices(images, angle, xbox, ybox):

    mat_dim = images * 2 + 1
    x_shift = np.zeros((mat_dim, mat_dim))
    y_shift = np.zeros((mat_dim, mat_dim))

    # shift in x and y direction due to angle
    x_comp = np.cos(angle*np.pi/180)*ybox
    y_comp = np.sin(angle*np.pi/180)*ybox

    for i in range(images + 1):
        x_shift[images, images + i] = i*xbox
        x_shift[images, images - i] = -i*xbox

    for i in range(1, images + 1):
        for j in range(images + 1):
            x_shift[images + i, images + j] = -i*x_comp + j*xbox
            x_shift[images - i, images + j] = i*x_comp + j*xbox
            x_shift[images + i, images - j] = -i*x_comp - j*xbox
            x_shift[images - i, images - j] = i*x_comp - j*xbox

    for i in range(1, images + 1):
        y_shift[images + i, :] = - i * y_comp
        y_shift[images - i, :] = i * y_comp

    return x_shift, y_shift


def duplicate_periodically(pts, box):
    """
    Duplicate points periodically in the +/- xyz directions. Assumes y box vector is angled with respect to x and
    the z vector points straight up
    :param pts: pts to be duplicated
    :param box: unitcell vectors in mdtraj format (use t.unitcell_vectors)
    :return: periodically extended system
    """

    nT = pts.shape[0]
    npts = pts.shape[1]
    p = np.zeros([nT, npts*27, 3])

    p[:, :npts, :] = pts

    for t in range(nT):

        p[t, npts:2*npts, :] = pts[t, :, :] + box[t, 2, :]
        p[t, 2*npts:3*npts, :] = pts[t, :, :] - box[t, 2, :]

        for i in range(3, 6):
            p[t, i*npts:(i+1)*npts, :] = p[t, (i-3)*npts:(i-2)*npts] + box[t, 0, :]

        for i in range(6, 9):
            p[t, i*npts:(i+1)*npts, :] = p[t, (i-6)*npts:(i-5)*npts] - box[t, 0, :]

        for i in range(9, 18):
            p[t, i*npts:(i+1)*npts, :] = p[t, (i-9)*npts:(i-8)*npts] + box[t, 1, :]

        for i in range(18, 27):
            p[t, i*npts:(i+1)*npts, :] = p[t, (i-18)*npts:(i-17)*npts] - box[t, 1, :]

    return p


def pbcs(pts, images, angle, xbox, ybox, frame):

    x_shift, y_shift = shift_matrices(images, angle, xbox, ybox)
    mat_dim = 2 * images + 1

    tot_pts = np.shape(pts)[1]

    translated_pts = np.zeros([3, mat_dim**2, tot_pts])

    if len(pts.shape) == 3:
        for p in range(tot_pts):
            for i in range(mat_dim):
                for j in range(mat_dim):
                    translated_pts[0, i*mat_dim + j, p] = x_shift[i, j] + pts[frame, p, 0]  # changed order for xlink.py and compatibility with mdtraj
                    translated_pts[1, i*mat_dim + j, p] = y_shift[i, j] + pts[frame, p, 1]
                    translated_pts[2, i*mat_dim + j, p] = pts[frame, p, 2]  # z position unchanged
    else:
        for p in range(tot_pts):
            for i in range(mat_dim):
                for j in range(mat_dim):
                    translated_pts[0, i*mat_dim + j, p] = x_shift[i, j] + pts[p, 0] / 10
                    translated_pts[1, i*mat_dim + j, p] = y_shift[i, j] + pts[p, 1] / 10
                    translated_pts[2, i*mat_dim + j, p] = pts[p, 2] / 10  # z position unchanged

    return translated_pts


if __name__ == "__main__":

    args = initialize()
    images = int(args.images)
    angle = float(args.angle)
    frame = 0

    gro = md.load(args.file)
    pos = gro.xyz
    full_box = gro.unitcell_vectors

    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                   full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]

    periodic_box = [3*float(i) for i in box_gromacs]

    id = [a.name for a in gro.topology.atoms]*27
    res = [a.residue.name for a in gro.topology.atoms]*27

    pt_periodic = duplicate_periodically(pos, full_box)
    pt_periodic[0, :, :] += [box_gromacs[0] + box_gromacs[5], box_gromacs[1], box_gromacs[2]]

    file_rw.write_gro_pos(pt_periodic[0, :, :], 'periodic.gro', box=periodic_box, ids=id, res=res)

    exit()

    pts = np.shape(pt_periodic)[2]
    duplicates = np.shape(pt_periodic)[1]

    all_positions = np.zeros([3, pts*duplicates])
    for i in range(duplicates):
        for j in range(pts):
            all_positions[:, i*pts + j] = pt_periodic[:, i, j]

    f = open('%s' % args.output, 'w+')

    f.write("This is a .gro file\n")
    f.write("%s\n" % (pts*duplicates))
    count = 0
    count1 = 1
    for j in range(duplicates):
        print(j)
        for i in range(pts):
            row = str(pt_periodic[:, j, i])

            # for full system
            if id[i] == 'NA':
                res = 'NA'
                count += 1
            else:
                res = 'HII'
                if id[i - 1] == 'NA':
                    count += 1
                if i != 0 and i % 137 == 0:
                    count += 1
            f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(count, '%s' % res, '%s' % id[i], count1, pt_periodic[0, j, i],
                                                                       pt_periodic[1, j, i], pt_periodic[2, j, i]) + "\n")

            # f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:9.3f}'.format(count, 'NA' , 'NA', count, pt_periodic[0, j, i],
            #                                                             pt_periodic[1, j, i], pt_periodic[2, j, i]) + "\n")
            # count += 1
            # count1 += 1

            count1 += 1
            if count1 == 100000:
                count1 = 0

    f.write('{:>3.6f}{:>3.6f}{:>3.6f}\n'.format(xbox*(2*images + 1), ybox*(2*images + 1), zbox))

    f.close()

    # f = open('NA_periodic_1.txt', 'w')
    # for j in range(duplicates):
    #     for i in range(pts):
    #         row = str(pt_periodic[:, j, i])
    #         row = row.replace("[", "")
    #         row = row.replace("]", "")
    #         f.write(row + "\n")
    # f.close()

