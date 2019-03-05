#!/usr/bin/env python

import argparse
import numpy as np
from llclib import file_rw
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Cut a cross section from a .gro file')  # allow input from user

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of input .gro file')
    parser.add_argument('-t', '--traj', type=str, help='Name of trajectory if you want to slice up a trajectory')
    parser.add_argument('-o', '--out', default='cross', type=str, help='Name of output .gro file')
    parser.add_argument('-p', '--plane', default='xy', help='Plane to slice through (xy, yx, xz, xz, xy, yz all work)')
    parser.add_argument('-s', '--slice', default=0.5, type=float, help='fraction along x y or z (depending on choice of '
                                                                       'plane) coordinate to slice through')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    if args.traj:
        t = md.load(args.traj, top=args.gro)
        print('Trajectory loaded')
    else:
        t = md.load(args.gro)
        print('GRO loaded')

    pos = t.xyz

    # box = t.unitcell_vectors[0, :, :]
    # vz = np.linalg.norm(box[2, :])
    # vx = np.linalg.norm(box[0, :])
    # vy = np.linalg.norm(box[1, :])
    #
    # # Cut at an angle
    # z = np.array([1, 0])*vz  # points w.r.t z which will be cut through
    # x = np.array([0.5, 1])*vx  # point w.r.t. x which will be cut through
    # m = (z[0] - z[1]) / (x[0] - x[1])  # slope
    # b = z[0] - m*x[0]  # intercept
    #
    # keep = []
    # for i in range(pos.shape[1]):
    #     if pos[0, i, 2] < m*pos[0, i, 0] + b:
    #         keep.append(i)

    d = 2  # default
    if args.plane == 'xy' or args.plane == 'yx':
        d = 2
    elif args.plane == 'xz' or args.plane == 'zx':
        d = 1
    elif args.plane == 'yz' or args.plane == 'zy':
        d = 0
    else:
        print('Please enter a valid plane to slice through')
        exit()

    dbox = np.linalg.norm(t.unitcell_vectors[0, d, :])

    keep = []
    for i in range(np.shape(pos)[1]):
        if pos[-1, i, d] <= dbox*args.slice:
            keep.append(i)

    print('Slicing and Dicing...')
    sliced = t.atom_slice(keep)
    print('Writing final configuration')
    if args.traj:
        traj = md.formats.XTCTrajectoryFile('%s.xtc' % args.out, mode='w', force_overwrite=True)  # creat mdtraj TRR trajectory object
        time = np.linspace(0, 1000, t.n_frames)  # arbitrary times. Times are required by mdtraj
        traj.write(sliced.xyz, time=time, box=t.unitcell_vectors)  # write the trajectory in .trr format
        file_rw.write_gro(sliced, '%s.gro' % args.out)
    else:
        file_rw.write_gro(sliced, '%s.gro' % args.out)
