#!/usr/bin/python

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
import argparse
import numpy as np
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Duplicate points periodically in the x-y directions')

    parser.add_argument('-i', '--images', default=5, help='Number of periodic images')
    parser.add_argument('-x', '--xbox', default=8.65745, help='Length of x box vector')
    parser.add_argument('-y', '--ybox', default=8.65745, help='Length of y box vector - NOTE: not y-component, but length of'
                                                        'entire vector')
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


def pbcs(pts, images, angle, xbox, ybox, frame):

    x_shift, y_shift = shift_matrices(images, angle, xbox, ybox)
    mat_dim = 2 * images + 1

    tot_pts = np.shape(pts)[1]
    translated_pts = np.zeros([3, mat_dim**2, tot_pts])

    for p in range(tot_pts):
        for i in range(mat_dim):
            for j in range(mat_dim):
                translated_pts[0, i*mat_dim + j, p] = x_shift[i, j] + pts[0, p, frame]
                translated_pts[1, i*mat_dim + j, p] = y_shift[i, j] + pts[1, p, frame]
                translated_pts[2, i*mat_dim + j, p] = pts[2, p, frame]  # z position unchanged

    return translated_pts


if __name__ == "__main__":

    args = initialize()
    images = int(args.images)
    angle = float(args.angle)
    xbox = float(args.xbox)
    ybox = float(args.ybox)
    frame = 0

    import Get_Positions
    pos = Get_Positions.get_positions('wiggle.gro', 'NA', 'HII', 'no')[0]
    print 'positions got'
    # For full system
    # pos = np.load('pos_array612ns')
    # id = np.load('identity_array612ns')

    pt_periodic = pbcs(pos, images, angle, xbox, ybox, frame)
    #
    f = open('pts_periodic', 'w')
    np.save(f, pt_periodic)
    f.close()

    #
    # pt_periodic = np.load('pts_periodic')
    # print np.shape(pt_periodic)
    # exit()

    f = open('NA_periodic.gro', 'w')

    pts = np.shape(pt_periodic)[2]
    duplicates = np.shape(pt_periodic)[1]

    all_positions = np.zeros([3, pts*duplicates])
    for i in range(duplicates):
        for j in range(pts):
            all_positions[:, i*pts + j] = pt_periodic[:, i, j]

    f = open('NA_positions_5_images', 'w')
    np.save(f, all_positions)
    f.close()
    exit()
    f.write("This is a .gro file\n")
    f.write("%s\n" % (pts*duplicates))
    count = 1
    for j in range(duplicates):
       for i in range(pts):

           row = str(pt_periodic[:, j, i])

           # for full system
           # if id[0, i] == 'NA':
           #     res = 'NA'
           # else:
           #     res = 'HII'
           #
           # f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:9.3f}'.format(count, '%s' % res, '%s' % id[0, i], count, pt_periodic[0, j, i],
           #                                                             pt_periodic[1, j, i], pt_periodic[2, j, i]) + "\n")

           f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:9.3f}'.format(count, 'NA' , 'NA', count, pt_periodic[0, j, i],
                                                                       pt_periodic[1, j, i], pt_periodic[2, j, i]) + "\n")

           count += 1

    f.write('   0.00000   0.00000  0.00000')

    f.close()

    f = open('NA_periodic_1.txt', 'w')
    for j in range(duplicates):
        for i in range(pts):
            row = str(pt_periodic[:, j, i])
            row = row.replace("[", "")
            row = row.replace("]", "")
            f.write(row + "\n")
    f.close()

