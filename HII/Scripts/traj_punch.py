#! /usr/bin/env python

import mdtraj as md
import argparse
import numpy as np
from llclib import file_rw
import os


def initialize():

    parser = argparse.ArgumentParser(description='Measure Membrane Thickness')

    parser.add_argument('-g', '--gro', type=str, default='wiggle.gro', help='Path to input file')
    parser.add_argument('-t', '--traj', type=str, default='traj_whole.xtc', help='Trajectory file (.xtc or .trr)')
    #parser.add_argument('-a', '--atoms', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', '0', 'O1', 'O2', 'O3', 'O4', 'NA', 'C7', 'C8', 'C9', 'C10', 'C21', 'C22', 'C23', 'C24', 'C35', 'C36', 'C37', 'C38'])
    #parser.add_argument('-a', '--atoms', default=['O5', 'O6', 'O7', 'O8', 'O9', 'O10', 'C46', 'C47', 'C48', 'C18', 'C19', 'C20', 'C32', 'C33', 'C34'])
    # parser.add_argument('-a', '--atoms', default=['NA'])
    parser.add_argument('-a', '--atoms', nargs='+', default=['C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20',
                                                   'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34',
                                                   'C35', 'C36', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48'])
    #parser.add_argument('-a', '--atoms', nargs='+')
    parser.add_argument('-discard', nargs='+', help='Atoms to not keep')
    parser.add_argument('-o', '--output', type=str, default='out')

    args = parser.parse_args()

    return args

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def punch():

    with open('%s/punch.txt' % location) as f:
        for line in f:
            print(line, end='')


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)
    print("Trajectory Loaded")

    if args.discard:
        keep = [a.index for a in t.topology.atoms if a.name not in args.discard]
    else:
        keep = [a.index for a in t.topology.atoms if a.name in args.atoms]

    new = t.atom_slice(keep)
    print("Trajectory PUNCHED!!!!!!!!!")
    punch()
    print("Saving...", end="")
    new.save_xtc('%s.xtc' % args.output)
    file_rw.write_gro(new, '%s.gro' % args.output)
    print("Done!")