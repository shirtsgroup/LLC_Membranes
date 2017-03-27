#! /usr/bin/env python

import argparse
from Last_Frame import extract_last_gro
from Thickness import thickness
from Get_Positions import get_positions as gp
import numpy as np
from llclib import physical


def initialize():

    parser = argparse.ArgumentParser(description = 'Run Cylindricity script')

    parser.add_argument('-i', '--input', default='wiggle_traj.gro')
    parser.add_argument('-t', '--traj', default='wiggle.trr', help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Path to input file')
    parser.add_argument('-c', '--comp', default='NA', help='Name of component whose concentration is needed')
    parser.add_argument('-b', '--buffer', default=0.1, help='percent depth into membrane where measurements are taken')
    parser.add_argument('-l', '--LC_type', default='NAcarb11V', help='Type of liquid crystal being studied')
    parser.add_argument('-s', '--solv', default='no', help='Is the system solvated or not?')

    args = parser.parse_args()

    return args

na = 6.022*10**23  # Avagadro's number


def conc(file, comp, b, lc, solv):

    extract_last_gro(file)  # this outputs a file called last_frame.gro
    thick, z_max, z_min = thickness('last_frame.gro')
    buffer = thick*float(b)
    z_max -= buffer
    z_min += buffer
    # Calculate concentration (an average of all frames)
    pos, _, _, box = gp(file, comp, lc, solv)
    no_comp = np.shape(pos)[1]
    trj_points = np.shape(pos)[2]
    count = np.zeros([trj_points])
    box_vol = np.zeros([trj_points])
    cross = np.zeros([trj_points])
    for t in range(trj_points):
        x_dim = np.sqrt(box[0][t]**2 + box[3][t]**2 + box[4][t]**2)
        y_dim = np.sqrt(box[1][t]**2 + box[5][t]**2 + box[6][t]**2)
        cross[t] = x_dim*y_dim
        box_vol[t] = x_dim*y_dim*thick
        for c in range(no_comp):
            if z_max >= pos[2, c, t] >= z_min:
                count[t] += 1

    factor = 1/(1*10**-27)  # convert from ions/nm^3 to ions/m^3
    conc = np.zeros([trj_points])
    for c in range(trj_points):
        conc[c] = (count[c]/box_vol[c])*factor

    avg_conc = np.mean(conc)
    std = np.std(conc)
    avg_cross = np.mean(cross)
    return avg_conc, std, avg_cross, thick, z_max, z_min

if __name__ == '__main__':
    args = initialize()
    import time
    start = time.time()
    # Avg, std, cross, thick, z_max, z_min = conc('%s' % args.input, '%s' %args.comp, '%s' % args.buffer, '%s' % args.LC_type, '%s' % args.solv)
    Avg, std, cross, thick, z_max, z_min = physical.conc('%s' % args.traj, '%s' %args.gro, '%s' % args.comp,
                                                         '%s' % args.buffer, '%s' % args.LC_type, '%s' % args.solv)

    print Avg
    print std
    print cross
    print z_max
    print z_min
    print 'Average Concentration: %s +/- %s mol/m^3' % (Avg, std)
    print 'Done in %s seconds' % (time.time() - start)