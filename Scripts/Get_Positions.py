#!/usr/bin/python

"""
This script basically deconstructs and organizes a .gro file into usable lists

Inputs:
    (i) A single frame or trajectory in .gro format
    (ii) The component you wish to study
    (iii) The type of liquid crystal (e.g. HII or BCC). It should have its attributes defined in LC_class.py

Outputs:
    (i) Positions of component 'comp' for each trajectory frame
    (ii) Velocities of component 'comp' for each trajectory frame
    (iii) Box vectors for each trajectory frame
    (iv) The timestamps on each trajectory frame
"""

import argparse
import LC_class

parser = argparse.ArgumentParser(description = 'Run Cylindricity script')

parser.add_argument('-i', '--input', default='wiggle_traj.gro', help='Name of file containing coordinates')
parser.add_argument('-c', '--component', default='NA', help = 'Path to input file')
parser.add_argument('-l', '--LC_type', default='HII', help = 'Name of this LC system. Should match its name in'
                                                             'LC_class.py')
parser.add_argument('-t', '--traj', default='yes', help = 'Do you want to read a trajectory .gro file?')
parser.add_argument('-s', '--solv', default='yes', help = 'Is the system solvated?')

args = parser.parse_args()


def traj_time(line, a, box_vect):
    i = 0
    while a[line][i:(i + 2)] != 't=':
        i += 1
    for k in range(0, 3):
        box_vect[k].append(float(a[line - 1][k*10:(k + 1)*10]))
    if len(a[line - 1]) > 30:  # If the last six entries are zero, they are not included in the .gro file
        for k in range(3, 9):
            box_vect[k].append(float(a[line - 1][k*10:(k + 1)*10]))
    timestamp = (float((a[line][i + 2: len(a[line]) - 1])))
    return timestamp, box_vect


def get_positions(input_file, comp, lc, solv):

    f = open(input_file, 'r')

    a = []
    for line in f:
        a.append(line)

    pos = [[], [], []]  # atomic positions [x, y, z]
    vel = [[], [], []]  # atom velocities  [vx, vy, vz]
    traj_points = []  # The time staps of each trajectory frame
    box = [[], [], [], [], [], [], [], [], []]  # box vectors, see: http://manual.gromacs.org/current/online/gro.html

    exec 'residues = LC_class.%s.residues' % lc

    if comp not in residues:
        residues.append(comp)  # Just in case the component is not normally associated with the LC_class
    if solv == 'yes':
        residues.append('SOL')

    # find the first line containing coordinates
    coords = 0  # start at the top
    while str.strip(a[coords][5:10]) not in residues:
        if a[coords].count('trjconv') != 0:
            timestamp, box = traj_time(coords, a, box)
            traj_points.append(timestamp)
        coords += 1

    for i in range(coords, len(a) - 1):
        if str.strip(a[i][5:10]) in residues:
            if str.strip(a[i][10:15]) == comp:
                pos[0].append(float(a[i][20:28]))
                pos[1].append(float(a[i][28:36]))
                pos[2].append(float(a[i][36:44]))
                vel[0].append(float(a[i][44:52]))
                vel[1].append(float(a[i][52:60]))
                vel[2].append(float(a[i][60:68]))
        else:
            if a[i].count('trjconv') != 0:
                timestamp, box = traj_time(i, a, box)
                traj_points.append(timestamp)

    # Now let's clean up what will become the output

    # The final box vectors take up the first entry in the box list so move those entries to the end
    for i in range(0, len(box)):
        box[i].append(box[i][0])  # Put the first box entry for each vector at the end each list
        del box[i][0]  # delete the first entry

    # Put all the positions and velocities into separate lists for each frame
    no_res_per_frame = len(pos[0])/len(traj_points)
    pos_frame = []
    vel_frame = []
    for i in range(0, len(traj_points)):
        pos_frame.append([])
        vel_frame.append([])
        for k in range(0, 3):
            pos_frame[i].append([])
            vel_frame[i].append([])

    for i in range(0, len(traj_points)):
        for k in range(0, 3):
            for j in range(0, no_res_per_frame):
                pos_frame[i][k].append(pos[k][i*no_res_per_frame + j])
                vel_frame[i][k].append(vel[k][i*no_res_per_frame + j])

    return pos_frame, vel_frame, traj_points, box

if __name__ == '__main__':
    posi, veli, trj_times, box = get_positions('%s' % args.input, '%s' % args.component, '%s' % args.LC_type, '%s' % args.solv)
    print veli[1][0]