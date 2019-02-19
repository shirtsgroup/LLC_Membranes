#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from LLC_Membranes.llclib import physical, topology
import tqdm


def initialize():

    parser = argparse.ArgumentParser(description='Calculate mean squared displacement (MSD) and diffusion coefficient'
                                                 'for a specific residue or set of atoms in a trajectory')

    parser.add_argument('-t', '--trajectory', default='wiggle.trr', help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of .gro coordinate file')
    parser.add_argument('-r', '--residue', type=str, help='Name of residue whose diffusivity we want')
    parser.add_argument('-begin', default=0, type=int, help='First frame to read')
    parser.add_argument('-end', default=-1, type=int, help='Last frame to read')
    parser.add_argument('-skip', default=1, type=int, help='Skip every n frames')
    parser.add_argument('-b', '--build_monomer', default='NAcarb11V.gro', help='Name of annotated coordinate file '
                        'describing monomer structure.')
    parser.add_argument('-atoms', nargs='+', help='Name of atoms whose collective diffusivity is desired')
    parser.add_argument('-a', '--axis', default='z', type=str, help='Which axis to compute msd along')
    parser.add_argument('-n', '--index', default=False, type=int, help='Solute number whose ztrace you want. If this is'
                                                                       'not specified, all of them will be shown')
    parser.add_argument('-pr', '--pore_radius', default=1.48, type=float, help='Max distance from pore center a water '
                        'molecule can exist in order to be counted as inside the pore')
    parser.add_argument('-cmax', '--colorbar_max', default=1.5, type=float, help='Maximum radial distance that gets'
                                                                                 'registered on the colorbar')

    return parser


class ZTrace(object):

    def __init__(self, traj, gro, residue, build_monomer, axis, begin=0, end=-1, skip=1, npores=4):

        self.t = md.load(traj, top=gro)[begin:end:skip]

        self.axis = ['x', 'y', 'z'].index(axis.lower())

        self.residue = topology.Residue(residue)
        self.monomer = topology.LC(build_monomer)

        res_ndx = [a.index for a in self.t.topology.atoms if a.residue.name == residue]
        residue_atom_names = [a.name for a in self.t.topology.atoms if a.residue.name == residue]
        masses = [self.residue.mass[x] for x in residue_atom_names[:self.residue.natoms]]

        self.nres = len(res_ndx) // len(masses)

        self.com = physical.center_of_mass(self.t.xyz[:, res_ndx, :], masses)

        self.spline = None
        self.radial_distance = np.zeros([self.t.n_frames, self.nres])

    def locate_pore_centers(self, npts_spline=10, save=True):

        pore_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms]

        self.spline = physical.trace_pores(self.t.xyz[:, pore_atoms, :], self.t.unitcell_vectors, npts_spline,
                                           save=save)[0]

    def radial_distances(self, npores=4):

        for t in tqdm.tqdm(range(self.t.n_frames), unit=' Frames'):
            d = np.zeros([npores, self.nres])
            for p in range(npores):

                d[p, :] = physical.radial_distance_spline(self.spline[t, p, ...], self.com[t, ...],
                                                          self.t.unitcell_vectors[t, ...])

            self.radial_distance[t, :] = d[np.argmin(d, axis=0), np.arange(self.nres)]

    def plot_trace(self, nr, colormap='plasma', cmax=None):

        fig, ax = plt.subplots()

        rd = self.radial_distance[:, nr]
        trace = self.com[:, nr, self.axis]

        if cmax is None:
            cmax = np.amax(rd)

        points = np.array([self.t.time / 1000, trace]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Create a continuous norm to map from data points to colors
        norm = plt.Normalize(0, cmax)
        lc = LineCollection(segments, cmap=colormap, norm=norm)
        # Set the values used for colormapping
        lc.set_array(rd)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)  # plot

        # Set bounds on plot
        ax.set_ylim(np.amin(trace), np.amax(trace))
        ax.set_xlim(0, self.t.time[-1] / 1000)

        # labels + colorbar
        ax.set_ylabel('$z$-coordinate (nm)', fontsize=14)
        ax.set_xlabel('Time (ns)', fontsize=14)
        fig.colorbar(line, ax=ax)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    trace = ZTrace(args.trajectory, args.gro, args.residue, args.build_monomer, args.axis, begin=args.begin,
                   end=args.end, skip=args.skip)

    trace.locate_pore_centers(save=True)

    trace.radial_distances()

    if args.index:
        trace.plot_trace(args.index, cmax=args.colorbar_max)
    else:
        for i in range(trace.com.shape[1]):
            trace.plot_trace(i, cmax=args.colorbar_max)
