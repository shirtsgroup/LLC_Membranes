#!/usr/bin/env python

"""
Extract a frame from a GROMACS trajectory
"""

import argparse
import mdtraj as md
from LLC_Membranes.llclib import file_rw, physical


def initialize():

    parser = argparse.ArgumentParser(description='Extract a frame from a GROMACS trajectory')

    parser.add_argument('-t', '--trajectory', default='PR_nojump.xtc', help='Path to input file.')
    parser.add_argument('-g', '--gro', default='em.gro', help='Name of .gro coordinate file.')
    parser.add_argument('-f', '--frame', default=0, type=int, help='Frame number (starting at 0)')
    parser.add_argument('-o', '--out', help='Name of output file. By default will name it frame_number.gro')
    parser.add_argument('--nowrap', action="store_true", help="Don't wrap coordinates so that they are in the box.")

    return parser


if __name__ == "__main__":

    args = initialize().parse_args()

    print('Loading trajectory...', flush=True, end='')
    t = md.load(args.trajectory, top=args.gro)
    print('Done!')

    if args.frame < 0:
        frame = t.n_frames + args.frame - 1
    else:
        frame = args.frame

    if args.out is None:
        out = 'frame_%s.gro' % frame
    else:
        out = args.out

    names = [a.name for a in t.topology.atoms]
    res = [a.residue.name for a in t.topology.atoms]

    if not args.nowrap:
        coords = physical.wrap_box(t.xyz[frame, ...], t.unitcell_vectors[frame, ...])
    else:
        coords = t.xyz[frame, ...]

    file_rw.write_gro_pos(coords, out, ids=names, res=res, ucell=t.unitcell_vectors[frame, ...])

    print('%s written' % out)

