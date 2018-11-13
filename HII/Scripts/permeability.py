#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import top
from llclib import physical
from scipy import spatial
import tqdm
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate permeability coefficient using the collective diffusion'
                                                 'approach')
    parser.add_argument('-t', '--trajectory', default='PR.xtc', type=str, help='Trajectory')
    parser.add_argument('-g', '--gro', default='berendsen.gro', type=str, help='Name of .gro coordinate file')
    parser.add_argument('-buf', '--buffer', default=0.05, type=float, help='Defines z limits on membrane channel. If '
                        'the value of buffer is 0.1, then 0.1*channel length is taken off the top and bottom of the'
                        'z box dimension')
    parser.add_argument('-r', '--residue', default='SOL', help='Name of residue whose permeability coefficient we want')
    parser.add_argument('--itp', default='/usr/local/gromacs/share/gromacs/top/amber99.ff/tip3p.itp', help='Name of itp'
                        'describing topology of residue')
    parser.add_argument('-b', '--nboot', default=200, help='Number of bootstrap trials to be run')
    parser.add_argument('-f', '--frontfrac', default=0.05, type=float, help='Where to start fitting line on msd curve')
    parser.add_argument('-F', '--fracshow', default=.2, type=float, help='Percent of graph to show, also where to stop '
                        'fitting line during diffusivity calculation')
    parser.add_argument('-ref', default='NA', type=str, help='Name of atom(s) to use as reference for locating pore '
                                                             'centers')
    parser.add_argument('-radius', default=1, type=float, help='Radius from pore center to look for solute molecules')

    args = parser.parse_args()

    return args


def dn(z, zmax, zmin):
    """
    Calculate the displacement of all positions from frame to frame
    :param z: trajectory of positions [nframes, npoints, 1]  --> z-direction only
    :return: a trajectory n(t) describing the delta n at each time step
    """

    n = np.zeros([z.shape[0]])
    for t in range(1, n.shape[0]):
        for i in range(z.shape[1]):
            current = z[t, i]
            previous = z[t - 1, i]
            if zmax >= current >= zmin:
                if zmax >= previous >= zmin:
                    n[t] += current - previous
                elif previous > zmax:
                    n[t] += current - zmax
                elif previous < zmin:
                    n[t] += current - zmin
            elif current > zmax:
                if zmax >= previous >= zmin:
                    n[t] += zmax - previous
            elif current < zmin:
                if zmax >= previous >= zmin:
                    n[t] += zmin - previous

    n /= (zmax - zmin)  # divide by channel length

    return n


def dn2(z, zmax, zmin):

    n = np.zeros([z.shape[0]])
    for t in tqdm.tqdm(range(1, z.shape[0])):
        for i in range(z.shape[1]):  # find out how far each ion has moved in the z direction this frame
            curr = z[t, i]  # current z position
            prev = z[t - 1, i]  # previous z position
            if zmin <= curr <= zmax:  # check that the ion is in the channel over which we are measuring
                current_frame = curr  # If it, that is its position at the current frame
            elif curr > zmax:  # If the particle is above zmax, we are still potentially interested in it
                current_frame = zmax
            elif curr < zmin:  # If the particle is below zmin, we are still potentially interested in it
                current_frame = zmin
            if zmin <= prev <= zmax:  # Do the same checks on the particle from the previous frame
                prev_frame = prev
            elif prev > zmax:
                prev_frame = zmax
            elif prev < zmin:
                prev_frame = zmin
            n[t] += (current_frame - prev_frame)

    return n


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.trajectory, top=args.gro)

    nT = t.n_frames  # number of frames

    res = args.residue
    if res == 'SOL':  # mdtraj changes the name from SOL to HOH
        res = 'HOH'

    selection = [a.index for a in t.topology.atoms if a.residue.name == res]

    pos = t.xyz[:, selection, :]  # positions of all atoms of interest

    topology = top.Top(args.itp)  # read topology
    atoms_per_residue = topology.natoms  # number atoms in a single residue
    matoms = topology.atom_masses  # mass of the atoms in residue
    mres = np.sum(matoms)  # total mass of residue

    com = np.zeros([nT, pos.shape[1]//atoms_per_residue, 3])  # track the center of mass of each residue

    for f in range(nT):
        for i in range(com.shape[1]):
            w = (pos[f, i*atoms_per_residue:(i+1)*atoms_per_residue, :].T * matoms).T  # weight each atom in the residue by its mass
            com[f, i, :] = np.sum(w, axis=0) / mres  # sum the coordinates and divide by the mass of the residue

    # locate pore centers

    locater = [a.index for a in t.topology.atoms if a.name == args.ref]

    centers = physical.avg_pore_loc(4, t.xyz[-1, locater, :])  # only looking at the last frame

    tree = spatial.cKDTree(centers.T)

    keep = []
    for i in range(com.shape[1]):
        # see how far the centers of mass of each solute is from it's nearest neighbor pore center
        if tree.query(com[-1, i, :2])[0] <= args.radius:
            keep.append(i)

    # Define the pore length
    L = t.unitcell_vectors[-1, 2, 2]  # use the last frame
    zmax = L - (args.buffer*L)
    zmin = args.buffer*L

    nt = dn(com[:, keep, 2], zmax, zmin)




