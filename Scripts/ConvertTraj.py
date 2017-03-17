#! /usr/bin/env python

"""
Convert a trajectory file in .gro format into simple ascii format with plain atom names and their coordinates
"""
import argparse
import Get_Positions
import numpy as np
import mdtraj as md


def initialize():

    parser = argparse.ArgumentParser(description='Duplicate points periodically in the x-y directions')

    parser.add_argument('-f', '--file', default='wiggle.trr', help='trajectory file (.xtc or .trr)')
    parser.add_argument('-c', '--coord', default='wiggle.gro', help='Coordinate file')
    parser.add_argument('-a', '--atoms', default='sys', help='Name of atom(s) you want in the output. Put sys if you '
                                                             'want all of them')
    parser.add_argument('-o', '--output', default='traj.txt', help='Name of output file')
    parser.add_argument('-r', '--res', default='HII', help='Name of residue wanted')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame number to begin reading')
    parser.add_argument('-e', '--end', type=int, help='specify if you want to end somewhere other than the last frame')
    parser.add_argument('--single_frame', help='Will create a copy of the current frame to make an artificial traj',
                        action="store_true")
    parser.add_argument('--avg_dims', help='Print the average dimensions Lx, Ly, and Lz', action="store_true")

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    if args.single_frame:
        t = md.load('%s' % args.coord)
        f = open('%s' % args.coord, 'r')
        a = f.readlines()
        f.close()
        box = a[-1].split()
    else:
        t = md.load('%s' % args.file, top='%s' % args.coord)
        box = t.unitcell_vectors  # get the unit cell lengths
    # t = md.load('%s' % args.file, top='%s' % args.coord)
    pos = t.xyz
    id = np.array([a.name for a in t.topology.atoms])
    nT = pos.shape[0]
    atoms = pos.shape[1]
    time = t.time
    print 'Trajectory read'

    # This part needs improvement in terms of generalizability
    for i in range(id.shape[0]):
        if id[i].count('C') != 0:
            id[i] = 'C'
        elif id[i].count('H') != 0:
            id[i] = 'H'
        elif id[i].count('O') != 0:
            id[i] = 'O'

    begin = args.begin
    if args.end:
        end = args.end
    else:
        end = nT

    print 'Writing file ...'
    with open('%s' % args.output, 'w') as f:
        f.write('%s %s\n' % (atoms, (end - begin + 1)))
        if args.single_frame:
            for i in range(2):
                f.write('Frame %s\n' % (i + 1))
                f.write('{:^6}{:^8}{:^8}{:^8}'.format('Atom', 'x', 'y', 'z') + '\n')
                for j in range(atoms):
                    if id[j] != 'PI':
                        f.write('{:6s}{:8.3f}{:8.3f}{:8.3f}\n'.format(id[j], pos[0, j, 0], pos[0, j, 1], pos[0, j, 2]))
                f.write('{:10s}{:10s}{:10s}{:10s}{:10s}{:10s}{:10s}{:10s}{:10s}\n'.format(box[0], box[1], box[2]
                                                                                    ,box[3], box[4], box[5]
                                                                                    ,box[6], box[7], box[8]))
        else:
            for i in range(begin - 1, end):
                print 'Frame %s written' % i
                # f.write('Frame %s, time = %s ps\n' % ((i + 1), time[i]))
                f.write('Frame %s\n' % (i + 1))
                f.write('{:^6}{:^8}{:^8}{:^8}'.format('Atom', 'x', 'y', 'z') + '\n')
                for j in range(atoms):
                    if id[j] != 'PI':
                        f.write('{:6s}{:8.3f}{:8.3f}{:8.3f}\n'.format(id[j], pos[i, j, 0], pos[i, j, 1], pos[i, j, 2]))
                f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}\n'.format(box[i, 0, 0], box[i, 1, 1], box[i, 2, 2]
                                                                                    ,box[i, 0, 1], box[i, 2, 0], box[i, 1, 0]
                                                                                    ,box[i, 0, 2], box[i, 1, 2], box[i, 2, 0]))

    print 'Output file %s written' % args.output

    f.close()
    if args.avg_dims:
        if args.single_frame:
            xyz = box[:3]
        else:
            xyz = [np.mean(box[:, 0, 0]), np.mean(box[:, 1, 1]), np.mean(box[:, 2, 2])]
        with open('dims.txt', 'w') as f:
            f.write('{:10f}{:10f}{:10f}'.format(xyz[0], xyz[1], xyz[2]))

    print 'Average box dimensions written'
