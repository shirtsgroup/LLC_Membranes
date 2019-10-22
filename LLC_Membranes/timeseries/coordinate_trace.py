#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from LLC_Membranes.llclib import physical, topology, file_rw
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
    parser.add_argument('-b', '--build_monomer', default='NAcarb11V', help='Name of monomer used to build unit '
                                                                           'cell.')
    parser.add_argument('-atoms', nargs='+', help='Name of atoms whose collective diffusivity is desired')
    parser.add_argument('-a', '--axis', default='z', type=str, help='Which axis to compute msd along')
    parser.add_argument('-n', '--index', default=False, nargs='+', help='Solute number whose ztrace you want. If this'
                                                                        'is not specified, all of them will be shown')
    parser.add_argument('-pr', '--pore_radius', default=1.48, type=float, help='Max distance from pore center a water '
                        'molecule can exist in order to be counted as inside the pore')
    parser.add_argument('-cmax', '--colorbar_max', default=1.5, type=float, help='Maximum radial distance that gets'
                                                                                 'registered on the colorbar')
    parser.add_argument('-load', '--load', default=False, help='Name of pickle file to load')
    parser.add_argument('-s', '--savename', default='trace.pl', help='Name to save ZTrace object under.')

    return parser


class CoordinateTrace(object):

    def __init__(self, traj, gro, residue, build_monomer, axis, begin=0, end=-1, skip=1):
        """ Initialize the calculation of a center of mass coordinate trace along a specified axis.

        :param traj: name of GROMACS trajectory (.xtc or .gro)
        :param gro: name of GROMACS coordinate file (.gro)
        :param residue: name of residue to track
        :param build_monomer: name of  monomer used to build LLC membrane
        :param axis: axis along which to trace (x, y or z)
        :param begin: first frame to track
        :param end: last frame to track
        :param skip: only record data every **skip** frames

        :type traj: str
        :type gro: str
        :type residue: str
        :type build_monomer: str
        :type axis: str
        :type begin: int
        :type end: int
        :type skip: int
        """

        self.t = md.load(traj, top=gro)[begin:end:skip]
        self.time = self.t.time

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

    def locate_pore_centers(self, npts_spline=10, save=True, savename='spline.pl'):
        """ Fit a spline through the centers of the pores based on the pore defining atoms
        (see :ref:`annotation-table-lc`) of the monomer used to construct the unit cell.

        :param npts_spline: number of points in each pore of the spline
        :param save: save the spline
        :param savename: name of spline

        :type npts_spline: int
        :type save: bool
        :type savename: str
        """

        pore_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms and
                      a.residue.name in self.monomer.residues]

        self.spline = physical.trace_pores(self.t.xyz[:, pore_atoms, :], self.t.unitcell_vectors, npts_spline,
                                           save=save, savename=savename)[0]

    def radial_distances(self):
        """ Calculate each residue of interest's distance from the pore center at each frame
        """

        npores = self.spline.shape[1]
        for t in tqdm.tqdm(range(self.t.n_frames), unit=' Frames'):
            d = np.zeros([npores, self.nres])
            for p in range(npores):

                d[p, :] = physical.radial_distance_spline(self.spline[t, p, ...], self.com[t, ...],
                                                          self.t.unitcell_vectors[t, ...])

            self.radial_distance[t, :] = d[np.argmin(d, axis=0), np.arange(self.nres)]

    def plot_mean_radial_distance(self, cutoff=None):

        if cutoff is not None:
            avg = []
            for t in range(self.t.n_frames):
                ndx = np.where(self.radial_distance[t, :] > cutoff)[0]
                avg.append(self.radial_distance[t, ndx].mean())
        else:
            avg = self.radial_distance.mean(axis=1)

        plt.plot(self.time / 1000, avg, lw=2)
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Mean distance from pore center (nm)', fontsize=14)
        plt.tick_params(labelsize=14)
        plt.tight_layout()
        plt.show()

    def plot_trace(self, nr, colormap='plasma_r', cmax=None, savename=None, show=True):
        """ Plot the coordinate trace of chosen residue centers of mass colored according to its radial distance from
        the closest pore center

        :param nr: index or indices of solutes to display. Indices are relative to the total number of residues,\
        starting from 0. For example, to look at the first and third residues pass [0, 2]
        :param colormap: name of colormap to use
        :param cmax: max value of color scale
        :param savename: name of plot
        :param show: show plot when finished plotting

        :type nr: int or list of int
        :type colormap: str
        :type cmax: float
        :type savename: str
        :type show: bool
        """

        fig, ax = plt.subplots()

        if type(nr) is not list:
            nr = [nr]

        # Create a continuous norm to map from data points to colors
        if cmax is None:
            cmax = np.amax(self.radial_distance[:, nr])  # maximum of all

        norm = plt.Normalize(0, cmax)

        for i in nr:
            print(i + 1)

            rd = self.radial_distance[:, i]
            trace = self.com[:, i, self.axis]

            points = np.array([self.time / 1000, trace]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)

            lc = LineCollection(segments, cmap=colormap, norm=norm)
            # Set the values used for colormapping
            lc.set_array(rd)
            lc.set_linewidth(2)
            line = ax.add_collection(lc)  # plot

        # Set reasonable bounds on plot
        ymin = np.amin(self.com[:, nr, self.axis])
        ymax = np.amax(self.com[:, nr, self.axis])
        span = ymax - ymin
        ymax += .05 * span
        ymin -= .05 * span

        ax.set_ylim(ymin, ymax)
        ax.set_xlim(0, self.time[-1] / 1000)

        # labels + colorbar
        ax.set_ylabel('$z$-coordinate (nm)', fontsize=14)
        ax.set_xlabel('Time (ns)', fontsize=14)
        cbar = fig.colorbar(line, ax=ax)
        cbar.set_label('Radial distance from pore center (nm)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        
        plt.tight_layout()

        if savename is not None:
            plt.savefig(savename)

        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.load:
        trace = file_rw.load_object(args.savename)
    else:
        trace = CoordinateTrace(args.trajectory, args.gro, args.residue, args.build_monomer, args.axis,
                                begin=args.begin, end=args.end, skip=args.skip)

        trace.locate_pore_centers(save=True)

        trace.radial_distances()

        file_rw.save_object(trace, '%s' % args.savename)

    # trace.plot_mean_radial_distance(args.pore_radius)
    # exit()

    if args.index:
        args.index = [int(i) for i in args.index]
        trace.plot_trace(args.index, cmax=args.colorbar_max)
    else:
        for i in range(trace.com.shape[1]):
            trace.plot_trace(i, cmax=args.colorbar_max)

