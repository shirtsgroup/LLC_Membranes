#! /usr/bin/env python

import file_rw
import mdtraj as md
import numpy as np


def thickness(filename, ref_atoms, *traj):
    """
    :param filename: name of .gro file
    :param ref_atoms: atoms which thickness will be based on
    :param traj: trajectory of positions
    :return: trajectory of thicknesses or single thickness based on max/min z coordinate of reference atoms
    """

    if traj:

        traj = np.asarray(traj)[0]  # optional arguments of the form *args need to be convert back to numpy arrays
        nT = traj.shape[0]  # number of trajectory points

        thick = np.zeros([nT])
        z_max = np.zeros([nT])
        z_min = np.zeros([nT])
        for t in range(nT):
            z_max_t = max(traj[t, :, 2])
            z_min_t = min(traj[t, :, 2])
            thick[t] = z_max_t - z_min_t
            z_max[t] = z_max_t
            z_min[t] = z_min_t

    else:
        f = open(filename, "r")  # .gro file whose positions of Na ions will be read

        a = []  # list to hold lines of file
        for line in f:
            a.append(line)

        line = 0
        while a[line].count('HII') == 0:
            line += 1

        z = []  # list to hold z positions of all atoms

        while a[line].count('HII') != 0:
            if str.strip(a[line][11:15]) in ref_atoms:
                z.append(float(a[line][36:44]))
            line += 1

        z_max = max(z)
        z_min = min(z)
        thick = z_max - z_min

    return thick, z_max, z_min


def conc(t, t_comp, b):
    """
    Calculate the concentration of the specified component
    :param t: mdtraj trajectory object for system being studied
    :param comp: component which you want the concentration of
    :param b: buffer. distance into membrane to go before starting calculation
    :return: concentration
    """

    box = t.unitcell_vectors
    last = t.slice(-1)
    file_rw.write_gro(last, 'last_frame.gro')
    thick, z_max, z_min = thickness('last_frame.gro')
    buffer = thick*b
    z_max -= buffer
    z_min += buffer
    thick = z_max - z_min

    # Calculate concentration (an average of all frames)
    pos = t_comp.xyz
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


def avg_pore_loc(npores, pos):
    """
    :param no_pores: the number of pores in the unit cell
    :param pos: the coordinates of the component(s) which you are using to locate the pore centers
                      (numpy array with dimensions: [no frames, no components, xyz coordinates, ] or just
                      [no components, xyz coordinates] for a single frame)
    :return: numpy array containing the x, y coordinates of the center of each pore at each frame
    """

    # Find the average location of the pores w.r.t. x and y

    if len(pos.shape) == 3:  # multiple frames

        nT = np.shape(pos)[0]
        comp_ppore = np.shape(pos)[1] / npores

        p_center = np.zeros([2, npores, nT])

        for i in range(nT):
            for j in range(npores):
                for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                    p_center[:, j, i] += pos[i, k, :2]
                p_center[:, j, i] /= comp_ppore  # take the average

    elif len(pos.shape) == 2:  # single frame

        comp_ppore = pos.shape[1] / npores
        p_center = np.zeros([2, npores])

        for j in range(npores):
            for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                p_center[:, j] += pos[:2, k]
            p_center[:, j] /= comp_ppore

    else:
        return 'Please use a position array with valid dimensions'
        exit()

    return p_center


def p2p(p_centers, distances):
    """
    :param p_centers: the x, y locations of the pore centers in the format return from avg_pore_loc()
    :param distances: the number of distinct distances between pores
    :return: all of the pore to pore distances
    """
    nT = np.shape(p_centers)[2]
    p2ps = np.zeros([distances, nT])  # distances in the order 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
    for i in range(nT):
        # So ugly ... sadness :(
        p2ps[0, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 1, i])
        p2ps[1, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 2, i])
        p2ps[2, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 3, i])
        p2ps[3, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 2, i])
        p2ps[4, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 3, i])
        p2ps[5, i] = np.linalg.norm(p_centers[:, 2, i] - p_centers[:, 3, i])

    return p2ps


def limits(pos, pcenters):
    """
    Estimate the pore 'radius' based on the position of some component and it's maximum deviation from the pore center
    :param: pos: the positions of all atoms included in making the estimate
    :param: pcenters: the x,y positions of the pore centers for each frame
    :return: an approximate pore radius. Beyond which, we have entered the alkane region
    """

    nT = pcenters.shape[2]
    npores = pcenters.shape[1]
    natoms= pos.shape[1]
    atom_ppore = natoms / npores

    deviation = np.zeros([nT, npores, atom_ppore])
    for f in range(nT):
        for i in range(atom_ppore):
            for j in range(npores):
                deviation[f, j, i] = np.linalg.norm(pos[f, j*atom_ppore + i, :2] - pcenters[:, j, f])  # left off here

    deviation = np.reshape(deviation, (nT, natoms))
    fr = np.zeros([nT])
    frstd = np.zeros([nT])

    for i in range(nT):
        fr[i] = np.mean(deviation[i, :])  # + np.std(deviation[i, :]) # maybe?
        frstd[i] = np.std(deviation[i, :])

    return fr, frstd