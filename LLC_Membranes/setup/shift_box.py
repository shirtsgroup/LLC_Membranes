#!/usr/bin/env python

import argparse
import mdtraj as md
from llclib import file_rw
import tqdm


def initialize():

    parser = argparse.ArgumentParser(description='Shift box location for easier trajectory processing')

    parser.add_argument('-t', '--traj', help='Name of GROMACS trajectory file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of GROMACS coordinate file')
    parser.add_argument('-x', '--xshift', default=0, type=float, help='Distance to shift in x-direction')
    parser.add_argument('-y', '--yshift', default=0, type=float, help='Distance to shift in y-direction')
    parser.add_argument('-z', '--zshift', default=0, type=float, help='Distance to shift in z-direction')
    parser.add_argument('-o', '--output', default='shifted', help='Name of output file (excluding extension)')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to start shifting at. Only frames after begin will be'
                                                         'shifted and written out.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    if args.traj:
        t = md.load(args.traj, top=args.gro)[args.begin:]
    else:
        t = md.load(args.gro)

    pos = t.xyz
    box = t.unitcell_vectors
    res = [a.residue.name for a in t.topology.atoms]
    ids = [a.name for a in t.topology.atoms]
    box_gromacs = [box[0, 0, 0], box[0, 1, 1], box[0, 2, 2], box[0, 0, 1], box[0, 2, 0],
                   box[0, 1, 0], box[0, 0, 2], box[0, 1, 2], box[0, 2, 0]]

    print('Shifting box with vector [%.2f, %.2f, %.2f]' % (args.xshift, args.yshift, args.zshift))
    for f in tqdm.tqdm(range(pos.shape[0]), unit='frame'):
        pos[f, :, :] += [args.xshift, args.yshift, args.zshift]

    if args.traj:
        traj = md.formats.XTCTrajectoryFile('%s.xtc' % args.output, mode='w', force_overwrite=True)  # creat mdtraj TRR trajectory object
        traj.write(pos, time=t.time, box=t.unitcell_vectors)  # write the trajectory in .trr format
    else:
        file_rw.write_gro_pos(pos[0, :, :], '%s.gro' % args.output, ids=ids, res=res, box=box_gromacs)