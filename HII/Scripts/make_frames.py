#!/usr/bin/env python

import numpy as np
import argparse
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Measure Membrane Thickness')

    parser.add_argument('-g', '--gro', type=str, default='wiggle.gro', help='Path to input file')
    parser.add_argument('-n', '--nframes', default=2, type=int, help='number of frames to make')
    parser.add_argument('-o', '--output', type=str, default='frames', help='Output trajectory')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)
    box = t.unitcell_vectors
    new_box = np.zeros([args.nframes, 3, 3])
    pos = np.zeros([args.nframes, t.n_atoms, 3])
    for i in range(args.nframes):
        pos[i, :, :] = t.xyz
        new_box[i, :, :] = box[0, :, :]

    t.xyz = pos
    t.unitcell_vectors = new_box
    t.time = np.linspace(0, 1, args.nframes)

    t.save_xtc('%s.xtc' % args.output)