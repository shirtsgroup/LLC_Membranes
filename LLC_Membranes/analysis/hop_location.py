#!/usr/bin/env python

import argparse
from LLC_Membranes.llclib import file_rw
from LLC_Membranes.timeseries.forecast_ctrw import System
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable


def initialize():

    parser = argparse.ArgumentParser(description='Calculate radial distribution of a monomer component')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='traj_whole.xtc', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-r', '--residue', help='Name of residue whose '
                        'radial distribution function we want to calculate')
    parser.add_argument('-l', '--load', default=False, type=str, help='Name of compressed pickel file to load')
    parser.add_argument('-cmax', '--set_cmap_max', default=False, type=float, help='Set the maximum on the colorbar')
    parser.add_argument('-cmap', '--colormap', default='plasma', type=str, help='Name of color map. See '
                        'https://matplotlib.org/examples/color/colormaps_reference.html for others.')
    parser.add_argument('-b', '--bins', default=25, type=int, help='Number of histogram bins')

    return parser


class HopLocation(System):

    def __init__(self, traj, gro, res):

        super().__init__(traj, gro, res)

    def plot_hop_locations(self, show=False, bins=25, colormap='viridis', cmap_max=None):

        r_digitized = np.digitize(self.r, np.linspace(0, np.max(self.r), bins)) - 1  # bin all distances
        heights = [len(np.where(r_digitized == i)[0]) for i in range(bins)]  # number of counts per distance bin

        hop_lengths = []  # make all hop lengths into one long list
        for i in self.hop_lengths:
            hop_lengths += i

        hop_lengths = np.array(hop_lengths)

        # Average hop length as a function of distance from pore center
        heat = np.zeros([bins])
        for i in range(bins):
            ndx = np.where(r_digitized == i)[0]
            if len(ndx) > 0:
                heat[i] = np.mean(np.abs(hop_lengths[np.where(r_digitized == i)[0]]))

        fig, ax = plt.subplots()

        # define color map limits
        if cmap_max is not None:
            cmax = cmap_max
        else:
            cmax = np.max(heat)

        print('Colorbar maximum: %.3f' % cmax)

        norm = plt.Normalize(0, cmax)
        cmap = plt.cm.get_cmap(colormap)
        colors = cmap(norm(heat))

        bar_locations = np.linspace(np.min(self.r), np.max(self.r), bins)
        bar_width = bar_locations[1] - bar_locations[0]

        for i in range(bins):
            ax.bar(bar_locations[i], heights[i], bar_width, color=colors[i])

        # set up colorbar
        sm = ScalarMappable(cmap=cmap, norm=plt.Normalize(0, cmax))
        sm.set_array([])

        cbar = plt.colorbar(sm)
        cbar.set_label('Average hop length ($\AA$)', rotation=90, fontsize=14)
        cbar.ax.tick_params(labelsize=14)

        # Make plot pretty
        ax.set_xlabel('Distance from pore center', fontsize=14)
        ax.set_ylabel('Number of hops', fontsize=14)
        plt.tick_params(axis='both', labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.load:
        hops = file_rw.load_object(args.load)

    else:
        hops = HopLocation(args.traj, args.gro, args.residue)
        hops.calculate_solute_partition(r=1.5, spline=True, membrane_residue='HII', write_tcl=False)
        hops.hops_and_dwells(penalty=0.25, nframes_dwell=10, locations=True)
        file_rw.save_object(hops, 'hop_locations.pl')

    if args.set_cmap_max:
        cmax = args.set_cmap_max
    else:
        cmax = None

    hops.plot_hop_locations(show=True, cmap_max=cmax, colormap=args.colormap, bins=args.bins)
