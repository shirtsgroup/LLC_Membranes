#!/usr/bin/python

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

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load('%s' % args.file, top='%s' % args.coord)
    pos = t.xyz
    id = np.array([a.name for a in t.topology.atoms])
    box = t.unitcell_vectors  # get the unit cell lengths
    nT = t.xyz.shape[0]
    atoms = t.xyz.shape[1]
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
        for i in range(begin - 1, end):
            # f.write('Frame %s, time = %s ps\n' % ((i + 1), time[i]))
            f.write('Frame %s\n' % (i + 1))
            f.write('{:^6}{:^8}{:^8}{:^8}'.format('Atom', 'x', 'y', 'z') + '\n')
            for j in range(atoms):
                if id[j] != 'PI':
                    f.write('{:6s}{:8.3f}{:8.3f}{:8.3f}'.format(id[j], pos[i, j, 0], pos[i, j, 1], pos[i, j, 2]) + '\n')
            f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}\n'.format(box[i, 0, 0], box[i, 1, 1], box[i, 2, 2]
                                                                                ,box[i, 0, 1], box[i, 1, 0], box[i, 2, 0]
                                                                                ,box[i, 0, 2], box[i, 1, 2], box[i, 2, 0]))

    print 'Output file %s written' % args.output

    f.close()