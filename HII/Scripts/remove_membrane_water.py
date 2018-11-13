#! /usr/bin/env python

import argparse
import mdtraj as md
from llclib import file_rw


def initialize():

    parser = argparse.ArgumentParser(description='Write .mdp files and topology files for various types of simulation')

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file to modify')
    parser.add_argument('-l', '--gap', default=3, type=float, help='Width of gap between membrane layers')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)
    full_box = t.unitcell_vectors
    z = full_box[0, 2, 2]
    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                   full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]

    dry_system_ndx = [a.index for a in t.topology.atoms if a.residue.name != 'HOH']
    water_ndx = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']
    res = [a.residue.name for a in t.topology.atoms]
    ids = [a.name for a in t.topology.atoms]

    pos = t.xyz[0, :, :]

    upper_limit = z - (args.gap / 2)
    lower_limit = (args.gap / 2)

    keep = []
    for i in water_ndx:
        if pos[i, 2] < lower_limit or pos[i, 2] > upper_limit:  # keep waters outside the membrane
            keep.append(i)
            keep.append(i + 1)
            keep.append(i + 2)

    new_system = pos[(dry_system_ndx + keep), :]
    all_res = res[:len(dry_system_ndx)] + ['SOL']*len(keep)
    all_names = ids[:len(dry_system_ndx)] + ['OW', 'HW1', 'Hw2']*(len(keep)//3)

    file_rw.write_gro_pos(new_system, 'dry_membrane.gro', ids=all_names, res=all_res, box=box_gromacs)