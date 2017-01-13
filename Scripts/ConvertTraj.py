#!/usr/bin/python

"""
Convert a trajectory file in .gro format into simple ascii format with plain atom names and their coordinates
"""
import argparse
import Get_Positions
import numpy as np


def initialize():

    parser = argparse.ArgumentParser(description='Duplicate points periodically in the x-y directions')

    parser.add_argument('-f', '--file', default='wiggle_traj.gro', help='File to replicate periodically')
    parser.add_argument('-a', '--atoms', default='sys', help='Name of atom(s) you want in the output. Put sys if you '
                                                             'want all of them')
    parser.add_argument('-o', '--output', default='traj.txt', help='Name of output file')
    parser.add_argument('-r', '--res', default='HII', help='Name of residue wanted')
    parser.add_argument('-F', '--frames', default='all', help='Number of frames to read (starting from beginning)')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    pos, _, _, box = Get_Positions.get_positions('%s' % args.file, '%s' % args.atoms, 'HII', 'no')
    print 'box and positions got'

    # Use the following block if you need to use the positions and/or box vectors again. It's especially useful if
    # Get_Positions takes a long time
    # f = open('box_array_thin', 'w')
    # np.save(f, box)
    # f.close()
    # f = open('pos_array_thin', 'w')
    # np.save(f, pos)
    # f.close()
    # box = np.load('box_array')
    # pos = np.load('pos_array612ns')

    no_comp = np.shape(pos)[1]

    f = open('%s' % args.file, 'r')
    a = []
    for line in f:
        a.append(line)
    f.close()

    id = np.zeros([1, no_comp], dtype=object)

    count = 2

    for j in range(no_comp):
        atom = str.strip(a[count][10:15])
        id[0, j] = atom
        count += 1

    print 'id array made'

    # again, save if you want
    # f = open('id_thin', 'w')
    # np.save(f, id)
    # f.close()
    # id = np.load('identity_array612ns')
    # id = np.load('id_plain')

    # This part needs improvement in terms of generalizability
    for i in range(np.shape(id)[1]):
        if id[0, i].count('C') != 0:
            id[0, i] = 'C'
        elif id[0, i].count('H') != 0:
            id[0, i] = 'H'
        elif id[0, i].count('O') != 0:
            id[0, i] = 'O'

    # again, save if you want
    # f = open('id_thin_plain', 'w')
    # np.save(f, id)
    # f.close()

    f = open('%s' % args.output, 'w')

    if args.frames == 'all':
        frames = np.shape(pos)[2]
    else:
        frames = int(args.frames)

    atoms = np.shape(id)[1]

    for i in range(frames):
        f.write('Frame %s' % (i + 1) + '\n')
        f.write('{:^6}{:^8}{:^8}{:^8}'.format('Atom', 'x', 'y', 'z') + '\n')
        for j in range(atoms):
            f.write('{:6s}{:8.3f}{:8.3f}{:8.3f}'.format(id[0, j], pos[0, j, i], pos[1, j, i], pos[2, j, i]) + '\n')
        f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}'.format(box[0, i], box[1, i], box[2, i], box[3, i]
                                                                                , box[4, i], box[5, i], box[6, i]
                                                                                , box[7, i], box[8, i]) + '\n')

    f.close()