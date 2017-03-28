#! /usr/bin/env python

import file_rw
import mdtraj as md
import numpy as np


def thickness(filename):

    f = open(filename, "r")  # .gro file whose positions of Na ions will be read

    a = []  # list to hold lines of file
    for line in f:
        a.append(line)

    line = 0
    while a[line].count('HII') == 0:
        line += 1

    z = []  # list to hold z positions of all atoms

    benz_carbs = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']

    while a[line].count('HII') != 0:
        if str.strip(a[line][11:15]) in benz_carbs:
            z.append(float(a[line][36:44]))
        line += 1

    z_max = max(z)
    z_min = min(z)
    thick = z_max - z_min

    return thick, z_max, z_min


def conc(traj, gro, comp, b, lc, solv):

    t = md.load('%s' % traj, top='%s' % gro)
    box = t.unitcell_vectors
    last = t.slice(-1)
    file_rw.write_gro(last, 'last_frame.gro')
    thick, z_max, z_min = thickness('last_frame.gro')
    buffer = thick*float(b)
    z_max -= buffer
    z_min += buffer
    thick = z_max - z_min

    # Calculate concentration (an average of all frames)
    keep = [a.index for a in t.topology.atoms if a.name == 'NA']
    comp_only = t.atom_slice(keep)
    pos = comp_only.xyz
    ncomp = pos.shape[1]  # number of components in the simulation which you want the concentration of
    nT = pos.shape[0]
    count = np.zeros([nT])
    box_vol = np.zeros([nT])
    cross = np.zeros([nT])
    for t in range(nT):
        x_dim = np.linalg.norm(box[t, 0, :])
        y_dim = np.linalg.norm(box[t, 1, :])
        cross[t] = x_dim*y_dim
        box_vol[t] = x_dim*y_dim*thick
        for c in range(ncomp):
            if z_max >= pos[t, c, 2] >= z_min:
                count[t] += 1

    factor = 1/(1*10**-27)  # convert from ions/nm^3 to ions/m^3. Trouble here for cython. Need to declare types
    conc = np.zeros([nT])
    for c in range(nT):
        conc[c] = (count[c]/box_vol[c])*factor

    avg_conc = np.mean(conc)
    std = np.std(conc)
    avg_cross = np.mean(cross)

    return avg_conc, std, avg_cross, thick, z_max, z_min