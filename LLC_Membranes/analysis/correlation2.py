#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from LLC_Membranes.analysis import detect_peaks
from LLC_Membranes.llclib import physical, topology
import tqdm
import matplotlib.pyplot as plt
from matplotlib import ticker
from LLC_Membranes.llclib import fast_rotate
from scipy.optimize import curve_fit
from scipy.interpolate import RegularGridInterpolator
import sys


def initialize():

    parser = argparse.ArgumentParser(description='Calculate and plot slice of full 3D correlation function')

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

    def __init__(self, gro, trajectory=None, atoms=None, res=None, com=True, bins=[100, 100, 100], begin=0, end=-1,
                 theta=120):
        """
        :param gro: GROMACS coordinate file
        :param trajectory: GROMACS trajectory file (.xtc or .trr)
        :param atoms: Names of atoms to be included in calculation. If there are multiple groups of atoms which you want
        to keep separate for the COM calculation, each group should be contained within a list within this list. Default
        is 'all' which will include all atoms in the system in the calculation.
        :param res: Names of residue that each atom group is associated with. If there are two groups, there should be
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
            if res is None:
                sys.exit("You must specify either a group of atoms or residue. Neither have been provided therefore "
                         "I can't do anything")
        if res is None:
            print('No residue names are attached to the atom groups provided, guessing at residue names instead. This '
                  'is probably a bad idea unless there is only one residue.')
            residues = [a.residue.name for a in self.t.topology.atoms]
            res = list(set(residues))
            print('Guessed residues: %s' % res)

        res = [topology.Residue(r) for r in res]

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
        self.correlation3d, self.edges = self.calculate_correlation_function()

        self.bin_centers = []
        for d in range(3):
            self.bin_centers.append(np.array([self.edges[d][i] + ((self.edges[d][i + 1] - self.edges[d][i]) / 2)
                                              for i in range(self.correlation3d.shape[d])]))

    def com(self):
        """
        Calculate center of mass of groups of atoms. Assumes groups are sequentially numbered.
        :return: trajectory of center of mass coordinates
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
        """
        rescale coordinates so that cell dimensions are constant over the simulation
        returns rescaled coordinates and average length
        """

        dims = np.linalg.norm(self.t.unitcell_vectors, axis=2)  # magnitude of unitcell vectors at each frame
        self.box = np.average(dims, axis=0)
        a = self.box / dims

        for f in range(self.positions.shape[0]):
            for i in range(3):
                self.positions[f, :, i] *= a[f, i]

    def wrap_coordinates(self):

        zv = [0.0, 0.0, 0.0]  # zero vector

        # put all atoms inside box - works for single frame and multiframe
        for it in range(self.positions.shape[0]):  # looped to save memory
            self.positions[it, ...] = np.where(self.positions[it, ...] < self.box, self.positions[it, ...],
                                               self.positions[it, ...] - self.box)  # get positions in periodic cell
            self.positions[it, ...] = np.where(self.positions[it, ...] > zv, self.positions[it, ...],
                                               self.positions[it, ...] + self.box)

    def calculate_correlation_function(self):

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

        # invert structure factor back to real space
        return np.fft.ifftn(sf), edges

    def make_slice(self, axis, radius=0, plot=True, show=False, fit=True):
        """
        Take a slice of the 3d correlation function along a specified axis
        :param axis: axis along which to take slice, int chosen from (0, 1, 2)
        :param radius: average all off-center slices within this radius of the center
        :return:
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

    def plot_slice(self, axis, show=False, fit=True, peak_locations=[], limits=None):
        """
        Plot slice generated by self.make_slice().
        :param axis: axis along which slice was taken
        :param show: show plot at the end
        :param fit: fit a decaying exponential function to the peaks
        :param peak_locations: locations of peaks to be fit
        :param limits: x and y limits of plot of form : ([x_lower, x_upper], [y_lower, y_upper]). If you want matplotlib
        to choose a specific axis limit for you, leave the list for that axis blank. If you want matplotlib to
        automatically choose both axes, don't specify anything for this option
        :return:
        """

        x = self.bin_centers[axis]
        plt.plot(self.bin_centers[axis][1:], self.slice, linewidth=2, label='Raw Data')

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
        plt.axes().tick_params(labelsize=14)
        plt.legend(loc=1, prop={'size': 16})
        plt.tight_layout()
        plt.tight_layout()

        if show:
            plt.show()

    def angle_average(self, ucell=None, NBR=80, rmax=-1, zbins=-1, zmax=-1, plot=True, show=False, save=False):

        # a work in progress

        ES = RegularGridInterpolator((self.bin_centers[0], self.bin_centers[1], self.bin_centers[2]), self.sf,
                                     bounds_error=False)

        THETA_BINS_PER_INV_ANG = 20.
        MIN_THETA_BINS = 10  # minimum allowed bins
        RBINS = NBR

        a1 = self.unit_cell[0]
        a2 = self.unit_cell[1]
        a3 = self.unit_cell[2]

        b1 = (np.cross(a2, a3)) / (np.dot(a1, np.cross(a2, a3)))
        b2 = (np.cross(a3, a1)) / (np.dot(a2, np.cross(a3, a1)))
        b3 = (np.cross(a1, a2)) / (np.dot(a3, np.cross(a1, a2)))

        b_inv = np.linalg.inv(np.vstack((b1, b2, b3)))

        if zbins == -1:
            ZBINS = self.freq_z.shape[0]  # 400
        else:
            ZBINS = zbins

        if zmax == -1:
            ZMAX = self.freq_z[-1]
        else:
            ZMAX = zmax

        XR = (self.freq_x[-1] - self.freq_x[0])
        YR = (self.freq_y[-1] - self.freq_y[0])

        if rmax == -1:
            Rmax = min(XR, YR) / 2.0
            Rmax *= 0.95
        else:
            Rmax = rmax

        rarr, rspace = np.linspace(0.0, Rmax, RBINS, retstep=True)
        zar = np.linspace(-ZMAX, ZMAX, ZBINS)

        oa = np.zeros((rarr.shape[0], zar.shape[0]))

        circ = 2. * np.pi * rarr  # circumference

        for ir in range(rarr.shape[0]):

            NTHETABINS = max(int(THETA_BINS_PER_INV_ANG * circ[ir]),
                             MIN_THETA_BINS)  # calculate number of bins at this r
            thetas = np.linspace(0.0, np.pi * 2.0, NTHETABINS, endpoint=False)  # generate theta array

            t, r, z = np.meshgrid(thetas, rarr[ir], zar)  # generate grid of cylindrical points

            xar = r * np.cos(t)  # set up x,y coords
            yar = r * np.sin(t)

            pts = np.vstack((xar.ravel(), yar.ravel(), z.ravel())).T  # reshape for interpolation

            if ucell is not None:
                pts = np.matmul(pts, b_inv)

            oa[ir, :] = np.average(ES(pts).reshape(r.shape), axis=1)  # store average values in final array

        mn = np.nanmin(oa)
        oa = np.where(np.isnan(oa), mn, oa)

        rad_avg = np.average(oa)
        oa /= rad_avg  # normalize

        # set up data for contourf plot by making it symmetrical
        self.angle_averaged = np.append(oa[::-1, :], oa[1:], axis=0)  # SF
        self.r_angle_averaged = np.append(-rarr[::-1], rarr[1:])  # R
        self.z_angle_averaged = np.append(z[:, 0, :], z[1:, 0, :], axis=0)[0]  # Z

        if plot:
            fig, ax = plt.subplots()
            MIN = np.amin(self.angle_averaged)
            MAX = np.amax(self.angle_averaged)*.25
            lvls = np.linspace(MIN, MAX, 200)
            heatmap = ax.contourf(self.r_angle_averaged, self.z_angle_averaged, self.angle_averaged.T, levels=lvls, cmap='jet',
                            extend='max')
            plt.colorbar(heatmap)
            plt.xlabel('$q_r (\AA^{-1})$')
            plt.ylabel('$q_z (\AA^{-1})$')

        if save:
            plt.savefig('rzplot.png')
        if show:
            plt.show()


def autocorrelation(x, largest_prime=500):
    """ FFT based autocorrelation function, which is faster than numpy.correlate. Efficiency is key in order to avoid
    headaches.
    :param x : multidimensional numpy array of y values of an equispaced timeseries. [npoints, ntrajectories]
    :param largest_prime : the largest prime factor of array length allowed. The smaller the faster. 1.6M points takes
    about 5 seconds with largest_prime=1000. Just be aware that you are losing data by truncating. But 5-6 data points
    isn't a big deal usually.
    """

    # x -= np.mean(x, axis=0)  # subtract the mean (alternate, but not equivalent way of calculating autocorrelation)

    # l = 2*x.shape[0] - 1
    #
    # while largest_prime_factor(l) >= largest_prime or l % 2 == 0:
    #     l -= 1
    #
    # x = x[:(l + 1) // 2, :]

    length = x.shape[0]*2 - 1

    fftx = np.fft.fft(x, n=length)
    ret = np.fft.ifft(fftx * np.conjugate(fftx))
    ret = np.fft.fftshift(ret)

    # divide each point in the autocorrelation function by the number of counts
    auto = ret[length // 2:].real
    # auto = np.mean(auto, axis=1)

    n = auto.shape[0]

    # return auto / (x.var()*np.arange(n, 0, -1))  # if subtracting mean, this is right thing to return

    return auto / np.arange(n, 0, -1)


def mc_inv(D, ucell):

    a1 = ucell[0]
    a2 = ucell[1]
    a3 = ucell[2]

    b1 = (np.cross(a2, a3))/(np.dot(a1, np.cross(a2, a3)))
    b2 = (np.cross(a3, a1))/(np.dot(a2, np.cross(a3, a1)))
    b3 = (np.cross(a1, a2))/(np.dot(a3, np.cross(a1, a2)))

    b_inv = np.linalg.inv(np.vstack((b1, b2, b3)))
    Dnew = np.zeros_like(D)

    X = D[..., 0]
    Y = D[..., 1]
    Z = D[..., 2]

    for ix in range(D.shape[0]):
        Dnew[ix, 0:3] += X[ix]*b_inv[0]

    for iy in range(D.shape[0]):
        Dnew[iy, 0:3] += Y[iy]*b_inv[1]

    for iz in range(D.shape[0]):
        Dnew[iz, 0:3] += Z[iz]*b_inv[2]

    return Dnew


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
        # Then calculate structure factor and invert to get correlation function.
        g = Correlation(args.gro, args.traj, atoms=args.atoms, res=args.res, com=args.center_of_mass, bins=bins,
                        begin=args.begin, end=args.end)

        if args.slice:

            #axes = {'x': 0, 'y': 1, 'z': 2, 'X': 0, 'Y': 1, 'Z': 2}
            axes = ['x', 'y', 'z']
            print(axes.index(args.slice.lower()))
            g.make_slice(axes.index(args.slice.lower()), radius=.225, plot=False)  # slice along z axis

            plt.imshow(g.correlation3d[0, :, :].real)
            plt.show()
            exit()

            g.plot_slice(axes[args.slice], show=True, fit=False,
                         limits=([0, g.bin_centers[axes[args.slice]][g.bin_centers[axes[args.slice]].size // 2]], []))

        exit()

        correlation = np.zeros(bins)

        x = np.array([1, 0, 0])
        xnorm = np.linalg.norm(x)

    if args.load:

        if len(args.slice) == 2:
            corr = np.load('correlation_%s%s.npz' %(args.slice[0], args.slice[1]))
        else:
            corr = np.load('correlation_%s.npz' %(args.slice[0]))

        correlation = corr['correlation']
        edges = corr['edges']
        frames = corr['frames']
        ncom = corr['ncom']

    else:
        frames = t.n_frames

        if args.slice == 'xy' or args.slice == 'yx':
            for frame in tqdm.tqdm(range(frames), unit='Frame'):
                for p in range(npores):
                    for l in range(args.layers):
                        for a in range(monomers_per_layer):

                            pt = p*com_per_pore + l*monomers_per_layer + a  # index of center of mass reference point
                            point = periodic_pts[frame, pt, :]  # coordinates of center of mass reference point
                            v = pore_spline[frame, p*args.layers + l, :] - point  # vector from point to pore center
                            rot = fast_rotate.rotate_vector(periodic_pts[frame, :, :], v, x)  # rotate all points by angle
                            trans = rot - rot[pt, :]  # translate all point so reference point is at the origin
                            H, edges = np.histogramdd(trans, bins=bins, range=hist_range)
                            correlation += H
                            # middle_x = bins[0] // 2
                            # middle_y = bins[1] // 2
                            # plt.plot(correlation[middle_x, middle_y, :])
                            # plt.show()

        else:
            # for frame in tqdm.tqdm(range(frames), unit='Frame'):
            #     for p in range(npores):
            #         for l in range(args.layers):
            #             for a in range(monomers_per_layer):
            #                 pt = p*com_per_pore + l*monomers_per_layer + a  # index of center of mass reference point
            #                 translated = periodic_pts[frame, :, :] - periodic_pts[frame, pt, :]  # make pt the origin
            #                 H, edges = np.histogramdd(translated, bins=bins, range=hist_range)
            #                 correlation += H
                            # middle_x = bins[0] // 2
                            # middle_y = bins[1] // 2
                            # #plt.plot(H[middle_x, middle_y, :])
                            # plt.plot(correlation[middle_x, middle_y, :])
                            # plt.show()

            # convert unit cell -- can be made more general
            theta = 2*np.pi / 3  # hard-coded = bad
            ucell = np.array([[1, 0, 0], [np.cos(theta), np.sin(theta), 0], [0, 0, 1]])

            # convert monoclinic cell to cubic cell
            print("transforming coordinates to monoclinic cell (theta={0:f} deg)".format(theta*180.0/np.pi))
            center_of_mass[..., 1] /= np.sin(theta)
            center_of_mass[..., 0] -= center_of_mass[..., 1]*np.cos(theta)

            L = np.linalg.norm(t.unitcell_vectors, axis=2)

            locations, L = rescale(center_of_mass, L)  # make unit cell constant size

            # redefine xyz values of grid
            box = t.unitcell_vectors
            x = np.linspace(0, L[0], int(bins[0]))
            y = np.linspace(0, L[1], int(bins[1]))  # for converted square box, y box length is same as x
            z = np.linspace(0, L[2], int(bins[2]))

            # redefine bins
            xbin = x[1] - x[0]
            ybin = y[1] - y[0]
            zbin = z[1] - z[0]

            zv = [0.0, 0.0, 0.0]  # zero vector

            # put all atoms inside box - works for single frame and multiframe
            for it in range(locations.shape[0]):  # looped to save memory
                locations[it, ...] = np.where(locations[it, ...] < L, locations[it, ...], locations[it, ...] - L)  # get positions in periodic cell
                locations[it, ...] = np.where(locations[it, ...] > zv, locations[it, ...], locations[it, ...] + L)

            # fourier transform loop
            fft = np.zeros([x.size - 1, y.size - 1, z.size - 1])
            for frame in tqdm.tqdm(range(frames), unit='Frames'):
                H, edges = np.histogramdd(locations[frame, ...], bins=(x, y, z))
                fft += np.abs(np.fft.fftn(H)) ** 2

            fft /= frames  # average of all frames

            # invert the fourier transform back to real space
            fft_inverse = np.fft.ifftn(fft)

            # only works for z slice currently
            centers_x = [edges[0][i] + ((edges[0][i + 1] - edges[0][i]) / 2) for i in range(fft_inverse.shape[0])]
            centers_y = [edges[1][i] + ((edges[1][i + 1] - edges[1][i]) / 2) for i in range(fft_inverse.shape[1])]
            centers_z = np.array([edges[2][i] + ((edges[2][i + 1] - edges[2][i]) / 2) for i in range(fft_inverse.shape[2])])

            # middle_x = fft_inverse.shape[0] // 2
            # middle_y = fft_inverse.shape[1] // 2

            # restricted zdf
            r = 5
            avg = []
            n = 0
            for i, ix in enumerate(centers_x):
                for j, jy in enumerate(centers_y):
                    if np.linalg.norm([ix - np.cos(theta)*jy, jy*np.sin(theta)]) < r:
                        avg.append([i, j])
                        n += 1

            zdf = np.zeros_like(centers_z)
            for i in range(len(avg)):
                zdf += fft_inverse[avg[i][0], avg[i][1], :].real

            zdf = zdf[1:]  # get rid of giant spike at zero
            zdf /= zdf.mean()

            zdf_full = np.mean(fft_inverse, axis=(0, 1))[1:]
            zdf_full /= zdf_full.mean()

            plt.plot(centers_z[1:], zdf_full, linewidth=2)
            # plt.plot(centers_z[1:], zdf_full, linewidth=2)
            plt.xlim(0, 4)
            # plt.ylim(0, 2)

            start = 13
            end = np.argmin(np.abs(np.array(centers_z) - 4))

            #Fit decaying exponential to peaks of oscillating correlation function
            peaks = detect_peaks.detect_peaks(zdf_full[start:end], mpd=12, show=False)  # adjust mpd if number of peaks comes out wrong
            if args.offset:
                peaks = peaks[1::2]  # every other peak starting at the second peak
            # peaks = [17, 34, 54, 74, 159] # full
            # peaks = [17, 35, 54, 163]

            # peaks = [9,  29,  51,  73,  101, 125, 143, 173]  # layered 300K ordered
            # peaks = [32, 78, 125, 165]  # offset 300K
            # peaks = [31, 77, 119, 162]  # offset 280K
            # peaks = [10, 30, 58, 80, 100, 118, 138, 157]  # layered 300K disordered
            #peaks = [33, 82, 115, 148]  # offset 300K disordered
            # peaks = [32, 82, 133]  # disorder offset
            print(peaks)
            # if len(peaks) > 4:
            #     peaks = peaks[:4]
            peaks = np.array(peaks)
            plt.scatter(centers_z[peaks + 1], zdf_full[peaks], marker='+', c='r', s=200, label='Peak locations')

            period = 0.438
            p = np.array([2, 10])  # initial guess at parameters
            bounds = ([0, 0], [np.inf, np.inf])
            solp, cov_x = curve_fit(exponential_decay, centers_z[peaks], zdf_full[peaks], p, bounds=bounds)

            # plt.plot(centers1[start:], 1 + solp[0]*np.exp(-centers1[start:]/solp[1]))

            plt.plot(centers_z[start:end], exponential_decay(np.array(centers_z[start:end]), solp[0], solp[1]), '--',
                     color='black', label='Least squares fit')
            print('Correlation length = %1.2f +/- %1.2f angstroms' % (10*solp[1], 10*np.sqrt(cov_x[1, 1])))

            p = [2.5, 0.4, np.pi, 0.5]  # amplitude, period, phase shift, correlation length. Pick values above what is expected
            bounds = ([0, 0.4, 0, 0.4], [np.inf, 0.6, np.inf, np.inf])  # bounds on fit parameters
            #
            # solp, pcov = curve_fit(sinusoidal_decay, centers_z[start:end], zdf_full[(start - 1):(end-1)], p, bounds=bounds)
            # print(solp)
            # # plot fit
            # plt.plot(centers_z[start:end], sinusoidal_decay(np.array(centers_z[start:end]), solp[0], solp[1], solp[2], solp[3]), '--',
            #          c='black', label='Least squares fit')
            #
            # # plt.plot(centers1[start:], 1 + solp[0]*np.exp(-centers1[start:]))
            #
            # print('Correlation length = %1.2f +/- %1.2f angstroms' % (10*solp[3], 10*np.sqrt(pcov[3, 3])))
            # print('Oscillation Period = %1.3f +/- %1.3f angstroms' % (10*solp[1], 10*np.sqrt(pcov[1, 1])))

            plt.xlabel('z (nm)')
            plt.ylabel('Frequency')
            plt.savefig('r1.png')
            plt.show()
            exit()
        if len(args.slice) > 1:
            np.savez_compressed('correlation_%s%s' % (args.slice[0], args.slice[1]), correlation=correlation, edges=edges, frames=frames, ncom=ncom)
        else:
            np.savez_compressed('correlation_%s' % (args.slice[0]), correlation=correlation, edges=edges, frames=frames, ncom=ncom)

    # Normalize 3D correlation function
    normalization = frames * ncom**2 / np.prod(correlation.shape)
    correlation /= normalization

    centers1 = [edges[dimensions[0]][i] + ((edges[dimensions[0]][i + 1] - edges[dimensions[0]][i])/2) for i in range(correlation.shape[dimensions[0]])]

    if len(dimensions) == 1:
        # only works for z slice currently
        centers_x = [edges[0][i] + ((edges[0][i + 1] - edges[0][i])/2) for i in range(correlation.shape[0])]
        centers_y = [edges[1][i] + ((edges[1][i + 1] - edges[1][i])/2) for i in range(correlation.shape[1])]
        middle_x = correlation.shape[0] // 2
        middle_y = correlation.shape[1] // 2
        r = 0.225
        r = 0.5

        avg = []
        n = 0
        for i, ix in enumerate(centers_x):
            for j, jy in enumerate(centers_y):
                if np.linalg.norm([ix, jy]) < r:
                    avg.append([i, j])
                    n += 1

        zdf = np.zeros_like(centers1)
        for i in range(len(avg)):
            zdf += correlation[avg[i][0], avg[i][1], :]

        shift = 0
        while centers1[shift] < 0.2:
            shift += 1

        plt.figure()
        ft = np.abs(np.fft.fft(zdf[shift:] - np.mean(zdf[shift:])))**2
        freq = np.fft.fftfreq(zdf[shift:].size, centers1[1] - centers1[0])
        ndx = np.argsort(freq)
        freq = freq[ndx]
        ft = ft[ndx]
        plt.plot(freq, ft)

        zdf /= np.mean(zdf)
        zdf = zdf[shift:]
        centers1 = np.array(centers1[shift:])
        np.savez_compressed('zdf.npz', zdf=zdf, centers=centers1)

        plt.figure()
        plt.plot(centers1, zdf, label='Raw data')

        if args.fit:
            start = 5  # discard data points up to start.
            # #Fit decaying exponential to peaks of oscillating correlation function
            # peaks = detect_peaks.detect_peaks(zdf, mpd=12, show=False)  # adjust mpd if number of peaks comes out wrong
            # if args.offset:
            #     peaks = peaks[1::2]  # every other peak starting at the second peak
            # peaks = [  9,  29,  51,  73,  101, 125, 143, 173]  # layered 300K ordered
            # # peaks = [32, 78, 125, 165]  # offset 300K
            # # peaks = [31, 77, 119, 162]  # offset 280K
            # # peaks = [10, 30, 58, 80, 100, 118, 138, 157]  # layered 300K disordered
            # #peaks = [33, 82, 115, 148]  # offset 300K disordered
            # # peaks = [32, 82, 133]  # disorder offset
            # print(peaks)
            # # if len(peaks) > 4:
            # #     peaks = peaks[:4]
            #
            # plt.scatter(centers1[peaks], zdf[peaks], marker='+', c='r', s=200, label='Peak locations')
            #
            # period = 0.438
            # p = np.array([2, 10])  # initial guess at parameters
            # bounds = ([0, 0], [np.inf, np.inf])
            # solp, cov_x = curve_fit(exponential_decay, centers1[peaks], zdf[peaks], p, bounds=bounds)
            #
            # # plt.plot(centers1[start:], 1 + solp[0]*np.exp(-centers1[start:]/solp[1]))
            #
            # plt.plot(centers1[start:], exponential_decay(np.array(centers1[start:]), solp[0], solp[1]), '--',
            #          color='black', label='Least squares fit')
            # print('Correlation length = %1.2f +/- %1.2f angstroms' % (10*solp[1], 10*np.sqrt(cov_x[1, 1])))

            # # fit decaying sinusoidal function to data
            p = [2.5, 0.4, np.pi, 0.5]  # amplitude, period, phase shift, correlation length. Pick values above what is expected
            bounds = ([0, 0.4, 0, 0.4], [np.inf, 0.6, np.inf, np.inf])  # bounds on fit parameters

            solp, pcov = curve_fit(sinusoidal_decay, centers1[start:], zdf[start:], p, bounds=bounds)
            print(solp)
            # plot fit
            plt.plot(centers1[start:], sinusoidal_decay(centers1[start:], solp[0], solp[1], solp[2], solp[3]), '--',
                     c='black', label='Least squares fit')

            # plt.plot(centers1[start:], 1 + solp[0]*np.exp(-centers1[start:]))

            print('Correlation length = %1.2f +/- %1.2f angstroms' % (10*solp[3], 10*np.sqrt(pcov[3, 3])))
            print('Oscillation Period = %1.3f +/- %1.3f angstroms' % (10*solp[1], 10*np.sqrt(pcov[1, 1])))

        plt.xlabel('Z distance separation (nm)', fontsize=14)
        plt.ylabel('Count', fontsize=14)
        plt.axes().tick_params(labelsize=14)
        plt.tight_layout()
        plt.legend(loc=1, prop={'size': 16})
        plt.ylim(0, 1.2*np.amax(zdf))
        plt.tight_layout()
        plt.savefig('z_correlation.png')

        # plt.figure()
        # ft = np.abs(np.fft.fft(zdf[start:] - np.mean(zdf[start:]))) ** 2
        # binsize = centers1[1] - centers1[0]
        # freq_x = np.fft.fftfreq(len(zdf[start:]), d=binsize)
        # ndx = np.argsort(freq_x)
        # freq_x = freq_x[ndx]
        # ft = ft[ndx]
        # plt.plot(freq_x, ft)
        plt.show()

    if len(dimensions) > 1:

        if not args.angle_average:

            centers2 = [edges[1][i] + ((edges[dimensions[1]][i + 1] - edges[dimensions[1]][i])/2) for i in range(bins[1])]
            ax = sum(range(0, 3)) - sum(dimensions)  # see which axis is missing out of [0 1 2]. E.g. if we are looking at the
            # yz cross sections, then dimensions will equal [1 2]. We want to sum cross sections along the x-axis which is
            # also the 0 axis. So ax will equal 0
            twoD = np.mean(correlation, axis=ax)

            # remove center region of 2D histogram since it is brightest and contains no useful information
            # first find the center bins (index)

            for x in range(len(centers1)):
                for y in range(len(centers2)):
                    if np.linalg.norm([centers1[x], centers2[y]]) < args.block_radius:
                        twoD[x, y] = 0

            fig = plt.figure()
            ax1 = fig.add_subplot(111)

            if args.plot_range:

                extent = [float(i) for i in args.plot_range]
                dim_1_start = 0
                while centers1[dim_1_start] < extent[0]:
                    dim_1_start += 1
                dim_1_end = 0
                while centers1[dim_1_end] < extent[1]:
                    dim_1_end += 1
                dim_2_start = 0
                while centers1[dim_2_start] < extent[2]:
                    dim_2_start += 1
                dim_2_end = 0
                while centers1[dim_2_end] < extent[3]:
                    dim_2_end += 1

                twoD = twoD[dim_1_start:dim_1_end, dim_2_start:dim_2_end]
            else:
                extent = [centers1[0], centers1[-1], centers2[0], centers2[-1]]

            heatmap = ax1.imshow(twoD.T, extent=extent, cmap='jet', interpolation='gaussian')
            cbar = plt.colorbar(heatmap)
            cbar.ax.tick_params(labelsize=14)
            plt.xlabel('%s (nm)' % args.slice[0], fontsize=14)
            plt.ylabel('%s (nm)' % args.slice[1], fontsize=14)
            plt.gcf().get_axes()[0].tick_params(labelsize=14)
            ax = plt.axes()
            ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
            ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
            ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.5))
            plt.tight_layout()

            if args.save:
                plt.savefig('%s.png' % args.save)
            if not args.noshow:
                plt.show()

        else:

            x_centers = np.array([edges[0][i] + ((edges[0][i + 1] - edges[0][i])/2) for i in range(len(edges[0]) - 1)])
            y_centers = np.array([edges[1][i] + ((edges[1][i + 1] - edges[1][i])/2) for i in range(len(edges[1]) - 1)])
            z_centers = np.array([edges[2][i] + ((edges[2][i + 1] - edges[2][i])/2) for i in range(len(edges[2]) - 1)])

            # Fourier transform of 3d correlation function with angle averaging
            # fft = np.abs(np.fft.fftn(correlation))**2
            #
            # xbin = x_centers[1] - x_centers[0]
            # ybin = y_centers[1] - y_centers[0]
            # zbin = z_centers[1] - z_centers[0]
            #
            # freq_x = np.fft.fftfreq(x_centers.size - 1, d=xbin)
            # ndx = np.argsort(freq_x)
            # freq_x = freq_x[ndx]
            #
            # freq_y = np.fft.fftfreq(y_centers.size - 1, d=ybin)
            # ndy = np.argsort(freq_y)
            # freq_y = freq_y[ndy]
            #
            # freq_z = np.fft.fftfreq(z_centers.size - 1, d=zbin)
            # ndz = np.argsort(freq_z)
            # freq_z = freq_z[ndz]
            #
            # fft = fft[ndx, :, :]
            # fft = fft[:, ndy, :]
            # fft = fft[:, :, ndz]
            #
            # averaged, rfin, zfin = angle_average(freq_x, freq_y, freq_z, fft)
            #
            # # xlim = 0.4
            # # zlim = 0.4
            #
            # MAX = 0.001
            # MIN = 0
            #
            # NLEVELS = 200
            # lvls = np.linspace(MIN, MAX, NLEVELS)  # contour levels
            #
            # plt.figure()
            # cs = plt.contourf(rfin, zfin[0], averaged.T, levels=lvls, cmap='jet', extend='max')
            # plt.xlabel('$q_r$ ($\AA^{-1}$)')
            # plt.ylabel('$q_z$ ($\AA^{-1}$)')
            # # plt.gcf().get_axes()[0].set_ylim(-zlim, zlim)
            # # plt.gcf().get_axes()[0].set_xlim(-xlim, xlim)
            # plt.colorbar(format='%.1f')
            # plt.tight_layout()
            # plt.show()
            # exit()
            final, rfin, zfin = angle_average(x_centers, y_centers, z_centers, correlation)

            # convert to angstroms
            rfin *= 10
            zfin *= 10

            xlim = 15
            zlim = 15

            NLEVELS = 200

            lvls = np.linspace(np.amin(final), 0.27*np.amax(final), NLEVELS)  # contour levels
            #lvls = np.linspace(, 0.35*np.amax(final))
            plt.figure()
            cs = plt.contourf(rfin, zfin[0], final.T, levels=lvls, cmap='seismic', extend='max')
            plt.xlabel('$r$ ($\AA$)')
            plt.ylabel('$z$ ($\AA$)')
            plt.gcf().get_axes()[0].set_ylim(-zlim, zlim)
            plt.gcf().get_axes()[0].set_xlim(-xlim, xlim)
            plt.colorbar(format='%.1f')
            plt.tight_layout()
            plt.savefig('angle_averaged_correlation.png')

            if args.invert:

                # Fourier transform 2D angle averaged correlation function
                FT = np.abs(np.fft.fftn(final.T)) ** 2

                rbin = rfin[1] - rfin[0]
                freq_r = np.fft.fftfreq(rfin.size - 1, d=rbin)
                ndr = np.argsort(freq_r)
                freq_r = freq_r[ndr]

                zbin = zfin[0, 1] - zfin[0, 0]
                freq_z = np.fft.fftfreq(zfin[0].size - 1, d=zbin)
                ndz = np.argsort(freq_z)
                freq_z = freq_z[ndz]

                FT = FT[ndz, :]
                FT = FT[:, ndr]

                plt.figure()
                levels = np.linspace(0, 1000, 200)
                plt.contourf(freq_r, freq_z, FT, levels=levels, cmap='jet', extend='both')

            plt.show()
