#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from LLC_Membranes.analysis import detect_peaks
from LLC_Membranes.llclib import physical, topology
import tqdm
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator
import sys


def initialize():

    parser = argparse.ArgumentParser(description='Calculate and plot slice of full 3D correlation function. Currently,'
                                                 ' only slices in the z direction for a monoclinic unit cell are '
                                                 ' implemented.')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file. Make sure to '
                                                                            'preprocess with gmx trjconv -pbc whole')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-a', '--atoms', nargs='+', action='append', help='Name of atoms to calculate correlation '
                        'function with respect to. The center of mass will be used')
    parser.add_argument('-r', '--res', nargs='+', help='Residue to create correlation function with '
                        'respect to. Will use center of mass. Will override atoms')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Start frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='End frame')
    parser.add_argument('-bins', nargs='+', default=100, type=int, help='Integer or array of bin values. If more than'
                        'one value is used, order the inputs according to the order given in args.axis')
    parser.add_argument('-com', '--center_of_mass', action="store_true", help='Calculate based on center of mass of'
                                                                               'args.atoms')
    parser.add_argument('-s', '--slice', help='Slice to be visualized')

    parser.add_argument('--itp', default='/home/bcoscia/PycharmProjects/GitHub/HII/top/Monomer_Tops/NAcarb11V.itp')
    parser.add_argument('-m', '--monomers_per_layer', default=5, type=int, help='Number of monomers per layer')

    parser.add_argument('--layers', default=20, type=int, help='Number of layers')
    parser.add_argument('--range', nargs='+', help='Range for histogram plot. A list of the form:'
                        '[dimension 1 lower, dimension 1 upper, dimension 2 lower, dimension 2 upper ...]')
    parser.add_argument('--save', type=str, help='Follow flag with name to save plot by (Not including file extension)')
    parser.add_argument('--scale', type=float, default=1, help='Scale color bar by this factor')
    parser.add_argument('-br', '--block_radius', default=0.25, type=float, help='radius of circle around origin to zero out')
    parser.add_argument('--noshow', action="store_true", help='If this flag is entered, plots will not be shown')
    parser.add_argument('--load', action="store_true", help='load previously calculated correlation function')
    parser.add_argument('-pr', '--plot_range', nargs='+', help='range to plot. A list of the form: '
                        '[dimension 1 lower, dimension 1 upper, dimension 2 lower, dimension 2 upper ...]')

    parser.add_argument('-offset', action="store_true", help='If system is in offset configuration, correlation length'
                                                             'will be calculated using every other peak')
    parser.add_argument('-aa', '--angle_average', action="store_true", help='Angle average the 3D correlation function'
                                                                            'about z-axis')
    parser.add_argument('-invert', action="store_true", help='Invert correlation function at end to get structure'
                                                             'factor')
    parser.add_argument('-fit', action="store_true", help='Fit decaying exponential function to correlation function')

    return parser


def sinusoidal_decay(x, a, b, c, d):
    """
    :param p: parameters: [period of oscillations, correlation time, amplitude, phase shift]
    :param x: x axis values
    :return: exponential function that sinusoidially decays to one
    """

    return 1 - a*np.cos((2*np.pi/b)*x + c)*np.exp(-x/d)


def exponential_decay(x, a, L):

    return 1 + a*np.exp(-x/L)


class Correlation(object):

    def __init__(self, gro, trajectory=None, atoms=None, res=None, bins=[100, 100, 100], begin=0, end=-1, theta=120):
        """ Prepare trajectory for correlation function calculation

        :param gro: GROMACS coordinate file
        :param trajectory: GROMACS trajectory file (.xtc or .trr)
        :param atoms: Names of atoms to be included in calculation. If there are multiple groups of atoms which you want\
        to keep separate for the COM calculation, each group should be contained within a list within this list. Default\
        is 'all' which will include all atoms in the system in the calculation.
        :param res: Names of residue that each atom group is associated with. If there are two groups, there should be\
        two residue names, even if they are the same residue
        :param com: Calculate correlation function based on center of mass of groups
        :param bins: number of bins in each dimension for histogramming
        :param begin: First frame of coordinates to use
        :param end: last frame of coordinates to use
        :param theta: angle between xy box vectors

        :type gro: str
        :type trajectory: str
        :type atoms: str or list
        :type res: str or list
        :type com: bool
        :type bins: list or int
        :type begin: int
        :type end: int
        :type theta: float
        """

        # Initialize certain variables
        self.positions = None
        self.box = None
        self.bins = bins
        self.slice = None
        self.theta = theta
        self.correlation3d = None
        self.edges = None
        self.bin_centers = []

        # Load coordinates
        if trajectory is None:
            # if no trajectory is provided, do correlation function of single frame
            self.t = md.load(gro)
            print('Configuration loaded')
        else:
            # load selected frame of a trajectory (default is all the frames)
            self.t = md.load(trajectory, top=gro)[begin:end]
            print('Trajectory loaded')

        if atoms is None:
            if res is None:  # better to raise an exception
                sys.exit("You must specify either a group of atoms or residue. Neither have been provided therefore "
                         "I can't do anything")

        if type(atoms[0]) is not list:
            atoms = [atoms]

        if res is None:
            print('No residue names are attached to the atom groups provided, guessing at residue names instead. This '
                  'is probably a bad idea unless there is only one residue.')
            residues = [a.residue.name for a in self.t.topology.atoms]
            res = list(set(residues))
            print('Guessed residues: %s' % res)

        if type(res) is list:
            res = [topology.Residue(r) for r in res]
        else:
            res = [topology.Residue(res)]

        # get indices of atoms to be included in calculation, separated into groups if specified
        self.keep = []
        for i, grp in enumerate(atoms):
            self.keep.append([a.index for a in self.t.topology.atoms if a.name in grp and a.residue.name ==
                              res[i].name])

        # get mass of all atoms that will be included in calculation
        self.mass = []
        for i, grp in enumerate(atoms):
            self.mass.append([res[i].mass[x] for x in grp])

        # Calculate centers of mass
        self.com()

        if self.theta != 90:

            self.theta *= (np.pi / 180)
            # convert monoclinic cell to cubic cell
            print("Transforming coordinates from monoclinic to cubic unit cell")
            self.positions[..., 1] /= np.sin(self.theta)
            self.positions[..., 0] -= self.positions[..., 1]*np.cos(self.theta)

        self.rescale()  # scale coordinates so box dimensions are the same for each frame
        self.wrap_coordinates()  # put all atoms in the box (faster than physical.wrapbox since this is a cube)

    def com(self):
        """ Calculate center of mass of groups of atoms. Assumes groups are sequentially numbered.
        """

        ngrps = len(self.mass)  # number of separate groups of which to calculate coms's
        n = [len(self.keep[i]) // len(self.mass[i]) for i in range(ngrps)]  # number of atoms in each group
        ndx = [sum(n[:i]) for i in range(len(n))]

        self.positions = np.zeros([self.t.n_frames, sum(n), 3])

        for f in range(self.t.n_frames):
            for g in range(ngrps):
                    for i in range(n[g]):
                        start = i * len(self.mass[g])
                        end = (i + 1) * len(self.mass[g])
                        w = (self.t.xyz[f, self.keep[g][start:end], :].T * self.mass[g]).T
                        self.positions[f, ndx[g] + i, :] = np.sum(w, axis=0) / sum(self.mass[g])

    def rescale(self):
        """ rescale coordinates so that cell dimensions are constant over the simulation
        """

        dims = np.linalg.norm(self.t.unitcell_vectors, axis=2)  # magnitude of unitcell vectors at each frame
        self.box = np.average(dims, axis=0)
        a = self.box / dims

        for f in range(self.positions.shape[0]):
            for i in range(3):
                self.positions[f, :, i] *= a[f, i]

    def wrap_coordinates(self):
        """ Put all atoms inside of cubic unit cell
        """

        zv = [0.0, 0.0, 0.0]  # zero vector

        # put all atoms inside box - works for single frame and multiframe
        for it in range(self.positions.shape[0]):  # looped to save memory
            self.positions[it, ...] = np.where(self.positions[it, ...] < self.box, self.positions[it, ...],
                                               self.positions[it, ...] - self.box)  # get positions in periodic cell
            self.positions[it, ...] = np.where(self.positions[it, ...] > zv, self.positions[it, ...],
                                               self.positions[it, ...] + self.box)

    def calculate_correlation_function(self):
        """ Calculate 3D correlation function in two main steps:
        (1) Calculate 3D structure factor using a discrete Fourier transform
        (2) Invert the structure factor back to real space
        """

        # define bin edges
        x = np.linspace(0, self.box[0], int(self.bins[0]))
        y = np.linspace(0, self.box[1], int(self.bins[1]))  # for converted square box, y box length is same as x
        z = np.linspace(0, self.box[2], int(self.bins[2]))

        # calculate structure factor loop (Histograms then fourier transforms)
        sf = np.zeros([x.size - 1, y.size - 1, z.size - 1])
        for frame in tqdm.tqdm(range(self.positions.shape[0]), unit='Frames'):
            H, edges = np.histogramdd(self.positions[frame, ...], bins=(x, y, z))
            fft = np.fft.fftn(H)  # fourier transform grid after subtracting mean
            sf += (fft * fft.conjugate()).real

        sf /= self.positions.shape[0]  # average over all frames

        self.correlation3d = np.fft.ifftn(sf)  # invert structure factor back to real space
        self.edges = edges

        for d in range(3):
            self.bin_centers.append(np.array([self.edges[d][i] + ((self.edges[d][i + 1] - self.edges[d][i]) / 2)
                                              for i in range(self.correlation3d.shape[d])]))

    def make_slice(self, axis, radius=0, plot=False, show=False, fit=False):
        """ Take a slice of the 3d correlation function along a specified axis

        :param axis: axis along which to take slice, where x=0, y=1, z=2
        :param radius: average all off-center slices within this radius of the center

        :type axis: int
        :type radius: float
        """

        other_dimensions = [i for i in range(3) if i != axis]
        avg = []
        n = 0
        for i, ix in enumerate(self.bin_centers[other_dimensions[0]]):
            for j, jy in enumerate(self.bin_centers[other_dimensions[1]]):
                if np.linalg.norm([ix - np.cos(self.theta)*jy, jy*np.sin(self.theta)]) < radius:  # only works for z-slice
                    avg.append([i, j])
                    n += 1

        g = np.zeros_like(self.bin_centers[axis])
        for i in range(len(avg)):
            g += self.correlation3d[avg[i][0], avg[i][1], :].real

        g = g[1:]  # get rid of giant spike at zero

        self.slice = g / g.mean()
        if plot:
            self.plot_slice(axis, show=show, fit=fit)

    def plot_slice(self, axis, show=True, fit=False, peak_locations=[], limits=None):
        """ Plot slice generated by make_slice()

        :param axis: axis along which slice was taken
        :param show: show plot at the end
        :param fit: fit a decaying exponential function to the peaks
        :param peak_locations: locations of peaks to be fit
        :param limits: x and y limits of plot of form : ([x_lower, x_upper], [y_lower, y_upper]). If you want matplotlib\
        to choose a specific axis limit for you, leave the list for that axis blank. If you want matplotlib to\
        automatically choose both axes, don't specify anything for this option

        :type axis: int
        :type show: bool
        :type fit: bool
        :type peak_locations: list
        :type limits: list or tuple
        """

        x = self.bin_centers[axis]
        plt.plot(self.bin_centers[axis][1:], self.slice, linewidth=2)

        if limits:
            xlim = limits[0]
            ylim = limits[1]
            if xlim:
                plt.xlim(xlim[0], xlim[1])
            if ylim:
                plt.ylim(ylim[0], ylim[1])

        if fit:

            if not peak_locations:
                peaks = np.array(detect_peaks.detect_peaks(self.slice, mpd=12, show=False))
            else:
                peaks = np.array(peak_locations)

            if limits:
                peaks = [i for i in peaks if x[i] <= xlim[1]]

            print("Peak Indices: %s" % peaks)

            p = np.array([2, 10])  # initial guess at fit parameters [amplitude, correlation length]
            bounds = ([0, 0], [np.inf, np.inf])  # constrain parameters so they are positive
            solp, cov_x = curve_fit(exponential_decay, x[peaks], self.slice[peaks], p, bounds=bounds)

            print('Correlation length = %1.2f +/- %1.2f angstroms' % (10 * solp[1], 10 * np.sqrt(cov_x[1, 1])))

            plt.scatter(x[peaks], self.slice[peaks], marker='+', c='r', s=200, label='Peak locations')
            plt.plot(x, exponential_decay(x, solp[0], solp[1]), '--', color='black', label='Least square fit')

        plt.xlabel('Z distance separation (nm)', fontsize=14)
        plt.ylabel('Count', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.legend(loc=1, prop={'size': 16})
        plt.tight_layout()
        plt.tight_layout()

        if show:
            plt.show()

    # def angle_average(self, ucell=None, NBR=80, rmax=-1, zbins=-1, zmax=-1, plot=True, show=False, save=False):
    #
    #     # a work in progress
    #
    #     ES = RegularGridInterpolator((self.bin_centers[0], self.bin_centers[1], self.bin_centers[2]), self.sf,
    #                                  bounds_error=False)
    #
    #     THETA_BINS_PER_INV_ANG = 20.
    #     MIN_THETA_BINS = 10  # minimum allowed bins
    #     RBINS = NBR
    #
    #     a1 = self.unit_cell[0]
    #     a2 = self.unit_cell[1]
    #     a3 = self.unit_cell[2]
    #
    #     b1 = (np.cross(a2, a3)) / (np.dot(a1, np.cross(a2, a3)))
    #     b2 = (np.cross(a3, a1)) / (np.dot(a2, np.cross(a3, a1)))
    #     b3 = (np.cross(a1, a2)) / (np.dot(a3, np.cross(a1, a2)))
    #
    #     b_inv = np.linalg.inv(np.vstack((b1, b2, b3)))
    #
    #     if zbins == -1:
    #         ZBINS = self.freq_z.shape[0]  # 400
    #     else:
    #         ZBINS = zbins
    #
    #     if zmax == -1:
    #         ZMAX = self.freq_z[-1]
    #     else:
    #         ZMAX = zmax
    #
    #     XR = (self.freq_x[-1] - self.freq_x[0])
    #     YR = (self.freq_y[-1] - self.freq_y[0])
    #
    #     if rmax == -1:
    #         Rmax = min(XR, YR) / 2.0
    #         Rmax *= 0.95
    #     else:
    #         Rmax = rmax
    #
    #     rarr, rspace = np.linspace(0.0, Rmax, RBINS, retstep=True)
    #     zar = np.linspace(-ZMAX, ZMAX, ZBINS)
    #
    #     oa = np.zeros((rarr.shape[0], zar.shape[0]))
    #
    #     circ = 2. * np.pi * rarr  # circumference
    #
    #     for ir in range(rarr.shape[0]):
    #
    #         NTHETABINS = max(int(THETA_BINS_PER_INV_ANG * circ[ir]),
    #                          MIN_THETA_BINS)  # calculate number of bins at this r
    #         thetas = np.linspace(0.0, np.pi * 2.0, NTHETABINS, endpoint=False)  # generate theta array
    #
    #         t, r, z = np.meshgrid(thetas, rarr[ir], zar)  # generate grid of cylindrical points
    #
    #         xar = r * np.cos(t)  # set up x,y coords
    #         yar = r * np.sin(t)
    #
    #         pts = np.vstack((xar.ravel(), yar.ravel(), z.ravel())).T  # reshape for interpolation
    #
    #         if ucell is not None:
    #             pts = np.matmul(pts, b_inv)
    #
    #         oa[ir, :] = np.average(ES(pts).reshape(r.shape), axis=1)  # store average values in final array
    #
    #     mn = np.nanmin(oa)
    #     oa = np.where(np.isnan(oa), mn, oa)
    #
    #     rad_avg = np.average(oa)
    #     oa /= rad_avg  # normalize
    #
    #     # set up data for contourf plot by making it symmetrical
    #     self.angle_averaged = np.append(oa[::-1, :], oa[1:], axis=0)  # SF
    #     self.r_angle_averaged = np.append(-rarr[::-1], rarr[1:])  # R
    #     self.z_angle_averaged = np.append(z[:, 0, :], z[1:, 0, :], axis=0)[0]  # Z
    #
    #     if plot:
    #         fig, ax = plt.subplots()
    #         MIN = np.amin(self.angle_averaged)
    #         MAX = np.amax(self.angle_averaged)*.25
    #         lvls = np.linspace(MIN, MAX, 200)
    #         heatmap = ax.contourf(self.r_angle_averaged, self.z_angle_averaged, self.angle_averaged.T, levels=lvls, cmap='jet',
    #                         extend='max')
    #         plt.colorbar(heatmap)
    #         plt.xlabel('$q_r (\AA^{-1})$')
    #         plt.ylabel('$q_z (\AA^{-1})$')
    #
    #     if save:
    #         plt.savefig('rzplot.png')
    #     if show:
    #         plt.show()


if __name__ == "__main__":

    ndimensions = 3  # always do full 3d correlation function
    args = initialize().parse_args()

    if not args.load:

        if type(args.bins) is list:
            if len(args.bins) > 1:
                bins = np.array(args.bins)
            else:
                bins = np.array([args.bins[0]]*ndimensions)
        else:
            bins = np.array([args.bins]*ndimensions)

        if args.atoms is None and args.res is None:
            args.atoms = [['C', 'C1', 'C2', 'C3', 'C4', 'C5']]

        # Slice position array and calculate centers of mass of groups
        g = Correlation(args.gro, args.traj, atoms=args.atoms, res=args.res, bins=bins, begin=args.begin, end=args.end)

        # Then calculate structure factor and invert to get correlation function.
        g.calculate_correlation_function()

        if args.slice:

            axes = ['x', 'y', 'z']
            slice_ndx = axes.index(args.slice.lower())
            g.make_slice(slice_ndx, radius=.225, plot=False)  # slice along z axis

            g.plot_slice(slice_ndx, show=True, fit=False,
                         limits=([0, g.bin_centers[slice_ndx][g.bin_centers[slice_ndx].size // 2]], []))
