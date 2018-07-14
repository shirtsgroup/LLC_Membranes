#! /usr/bin/env python

""" Calculate structure factor of single 3D configurations or trajectories of configurations """

import tqdm
import argparse
import mdtraj as md
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import curve_fit


def initialize():

    parser = argparse.ArgumentParser(description='Create configurations in 1, 2 or 3 dimensions and take the discrete'
                                                 'fourier transform of the configuration')

    # Options that affect structure factor calculation
    parser.add_argument('-g', '--grid', nargs='+', default=[100, 100, 100], help='Number of real space grid points in'
                        'each direction (list of ints). Length of array must match number of dimensions')

    # Options for GROMACS trajectories (implementation needs to be copied over (and restructured) from fft3d.py)
    parser.add_argument('-gro', help='Name of coordinate file to fourier transform')
    parser.add_argument('-traj', help='Name of trajectory file to fourier transform')
    parser.add_argument('-avg', action="store_true", help='Valid when -traj is specified. Average histograms of atomic'
                        'positions over trajectory and then take fourier transform once of averaged array')
    parser.add_argument('-begin', default=0, type=int, help='Frame to begin calculations at')
    parser.add_argument('-end', default=-1, type=int, help='Frame to end calculations at')
    parser.add_argument('-noise', default=1, type=float, help='Frame to end calculations at')

    # Options for custom trajectories
    parser.add_argument('-box', nargs='+', default=[85, 85, 37], help='length of box vectors. Only orthorhombic boxes'
                                                                      'are implemented')
    parser.add_argument('-nframes', default=100, type=int, help='Number of frames to create')
    parser.add_argument('--random_columns', action="store_true", help='Create a trajectory with columns made of'
                        'equispaced points but are randomly displaced in the z-diretion with respect to other columns')
    parser.add_argument('--random_layers', action="store_true", help='Create hexagonally packed pores with layers '
                        'that are randomly rotated about the z-axis')
    parser.add_argument('-ncol', '--ncolumns', default=10, type=int, help='The number of columns in the x and y '
                        'dimensions. The total number of columns in the unit cell will be ncol^2. In the case of '
                        'hexagonal columns, this defines the number of pore centers.')
    parser.add_argument('-dbwl', default=3.7, type=float, help='Distance between layers (float, angstroms)')
    parser.add_argument('-nonoise', action="store_false", help='Turn off random column displacement')
    # hexagonally packed columns (combined with --random_columns flag)
    parser.add_argument('--hexagonal', action="store_true", help='Create a hexagonal packed columns that mimic the HII'
                                                                 'phase')
    parser.add_argument('--ncol_per_pore', default=5, type=int, help='Number of columns surrounding each pore center')
    parser.add_argument('--pore_radius', default=5, type=float, help='Distance to place each column from pore center')
    parser.add_argument('--cell_theta', default=120, type=float, help='Angle between vectors defining xy plane of'
                                                                      'monoclinic box')
    parser.add_argument('-npores', default=2, type=int, help='Number of pores in each dimension, similar to ncol')
    parser.add_argument('-thermal_disorder', default=[0, 0, 0], nargs='+', type=float, help='Degree of thermal noise in' 
                        ' each dimension expressed as a fraction of the distance between layers.')

    # The following are meant for custom trajectories but are not implemented. See fft3d.py for their implementation
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of layers in z direction (int)')
    parser.add_argument('-s', '--sigmas', default=[0.01, 0.01, 0.01], nargs='+', help='Sigma in each direction. (list'
                        'of floats). Length of list must match number of dimensions')
    parser.add_argument('-offset', action="store_true", help='Build offset configuration')
    parser.add_argument('-offset_angle', default='symmetric', help='Angle to offset adjacent layer (degrees). Default '
                        'is half of 360/args.nmon')
    parser.add_argument('-rd', '--radial_displaced', help='Distance to radially displaced layers')
    parser.add_argument('-rd_angle', help='Angle (degrees), w.r.t. x-axis, directing where layers will be radially'
                                          'displaced')
    parser.add_argument('-nmon', default=5, type=int, help='Number of monomers per layer (int)')
    parser.add_argument('-pore_radius', type=float, default=5, help='Pore radius (float)')
    parser.add_argument('-stagger', default=False, help='Deviation from uniform spacing of layers (float). e.g. If '
                        'layers are spaced 3.7 apart and you choose -stagger 1, then the first layer would be at zero,'
                        'the second at (3.7 - args.stagger) = 2.7, and the third at 2*3.7=7.4')
    parser.add_argument('-o', '--output', default='rzplot', help='Name of angle averaged plot to save. (no file extension)')
    parser.add_argument('--noshow', action="store_true", help='Do not show plots at end')

    return parser


def rotate_about_z(theta, xyz):

    origin = np.array([0, 0, 0])
    pos = np.copy(xyz)
    R = np.zeros([3, 3])
    R[0, :] = [np.cos(theta), -np.sin(theta), 0]
    R[1, :] = [np.sin(theta), np.cos(theta), 0]
    R[2, :] = [0, 0, 1]

    # translate points to origin
    center = np.mean(pos, axis=0)
    pos = translate(pos, center, origin)

    for i in range(np.shape(pos)[0]):
        pos[i, :] = np.dot(R, pos[i, :])

    return translate(pos, origin, center)


def translate(xyz, before, after):
    """
    :param xyz: coordinates of set of points to be translated [npts, 3]
    :param before: reference coordinate location before [3]
    :param after: reference coordinate location after [3]
    :return: translated points with respect to reference coordinate before/after locations [npts, 3]
    """

    pos = np.copy(xyz)
    direction = after - before

    translation = np.matrix([[1, 0, 0, direction[0]], [0, 1, 0, direction[1]],
                             [0, 0, 1, direction[2]], [0, 0, 0, 1]])

    b = np.ones([1])
    for i in range(pos.shape[0]):
        coord = np.concatenate((pos[i, :], b))
        x = np.dot(translation, coord)
        pos[i, :] = x[0, :3]

    return pos


def find_peaks(x, y, tol):
    """
    :param x: x values
    :param y: y values
    :param tol: determines what is considered a peak. The difference in height between a point and its neighbors must be
    at least "tol" times higher in order to be returned as a peak.
    :return:
    """

    peaks = []
    for i in range(1, x.size - 1):
        if y[i] > y[i - 1] and y[i] > y[i + 1]:
            if np.abs(y[i] - y[i - 1]) / y[i - 1] > tol and np.abs(y[i] - y[i + 1]) > tol:
                peaks.append(i)

    return np.array(peaks)


def lorentz(points, a, b, c):
    """
    :param p: lorentzian parameters : [full width half max (FWHM), position of maximum, maximum heigth]
    :param p: position
    :return:
    """

    w = a / 2

    x = (b - points) / w

    return (c / (np.pi*w)) / (1 + x ** 2)


def gaussian(points, mean, sigma, amplitude, yshift):

    return 1 + yshift + (amplitude / np.sqrt(2*np.pi*sigma**2)) * np.exp(-(points - mean)**2/(2*sigma**2))


def errorfunc(p, points, z):
    return lorentz(p, points) - z


def z_correlation(z, L, v=0.1):
    """
    Calculate where to place monomers on the z-axis so that a given correlation length is obtained
    :param z: mean z-positions where monomers will be placed with gaussian probability np.array([n_layers])
    :param L: desired correlation length [float]
    :param v: variance in z position of monomer head groups
    :return: locations [np.array[nlayers])
    """

    n = z.shape[0]
    cov = np.zeros([n, n])  # initialize covariance matrix

    decay = v*np.exp(-z / L)  # decay of covariance

    # decay[1:] += np.exp(-z[::-1][:-1]/L) # for periodicity (?)

    for i in range(z.shape[0]):
        cov[i, i:] = decay[:(n - i)]
        cov[i:, i] = decay[:(n - i)]

    # plt.imshow(cov, extent=[1, 20, 1, 20])
    # cbar = plt.colorbar()
    # cbar.set_ticks([0.02, 0.04, 0.06, 0.08, 0.1])
    # cbar.set_ticklabels([0.02, 0.04, 0.06, 0.08, 0.1])
    # ax = plt.gca()
    #
    # ticks = np.linspace(1, z.size - 1, z.size // 2, dtype=int)
    # ax.xaxis.set_ticks(ticks)
    # ax.yaxis.set_ticks(ticks)
    # plt.xlabel('Scatterer Number')
    # plt.ylabel('Scatterer Number')
    # plt.show()
    # exit()
    locations = np.random.multivariate_normal(z, cov)

    return locations


class Trajectory(object):

    def __init__(self):

        # initialize variables
        self.locations = None
        self.box = None
        self.nframes = 0
        self.sf = None
        self.freq_x = None
        self.freq_y = None
        self.freq_z = None
        self.slice = None
        self.unit_cell = None
        self.theta = 0
        self.atomic_form_factor = 0
        self.r_angle_averaged = 0
        self.z_angle_averaged = 0
        self.angle_averaged = 0

    def square_column_grid(self, ncolumns, npoints, frames=1, z_separation=3.7, xy_separation=1, bounds=None, noise=True):
        """
        :param ncolumns: Number of columns in 1 direction. There will be ncolumns**2 total columns
        :param npoints: Number of points in column array
        :param z_separation: distance between points in columns
        :param xy_separation: distance between columns in xy directions (same in both)
        :param bounds: bounds of histogram in each dimension
        :return: grid of locations
        """

        self.nframes = frames
        self.locations = np.zeros([self.nframes, ncolumns ** 2 * npoints, 3])

        # can use self.box if assumed that one corner of the box is at the origin
        if bounds:
            z_separation = bounds[2][1] / npoints
            xy_separation = bounds[0][1] / ncolumns  # assume x and y have equal dimensions (for now)

        x = np.linspace(0, (ncolumns - 1) * xy_separation, ncolumns)
        X, Y = np.meshgrid(x, x)
        column = np.linspace(0, z_separation * (npoints - 1), npoints)

        print('z-spacing: %.2f' % z_separation)
        print('xy-spacing: %.2f' % xy_separation)

        for t in range(self.nframes):
            for c in range(ncolumns ** 2):
                if noise:
                    shift = (z_separation / 2) * np.random.uniform(-1, 1)  # shift column vertically by a random amount
                else:
                    shift = 0
                self.locations[t, c * npoints:(c + 1) * npoints, 2] = column + shift
                self.locations[t, c * npoints:(c + 1) * npoints, :2] = [X[c % ncolumns, c // ncolumns],
                                                                Y[c % ncolumns, c // ncolumns]]

    def set_up_hexagonal(self, cell_theta):

        if self.box[0] != self.box[1]:
            print('WARNING: The x and y box lengths must be equal for this script to properly implement hexagonal '
                  'periodicity. Setting y box length equal to x box length.')
            self.box[1] = self.box[0]

        self.theta = cell_theta * np.pi / 180.0  # theta for monoclinic unit cell
        self.unit_cell = np.array([[1, 0, 0], [np.cos(self.theta), np.sin(self.theta), 0], [0, 0, 1]])

    def hexagonal_column_grid(self, npores, ncol_per_pore, r, npoints, frames=1, noise=True, thermal_disorder=[0, 0, 0]):

        self.nframes = frames
        self.locations = np.zeros([self.nframes, npores ** 2 * npoints * ncol_per_pore, 3])

        xy_pore_centers = np.zeros([npores**2, 2])
        dx = self.box[0] / npores  # distance between pores in x direction

        if r > (dx / 2):
            print('WARNING: The pore radius is such that pores intersect and  columns will be placed outside of the '
                  'unit cell which will disrupt periodicity. \nSetting r to %.2f. \nEither change the radius, change '
                  'the box dimensions, or change the number of pores in the unit cell.' % (dx/2))

        for i in range(npores):
            row_x = i*self.unit_cell[1, 0]*dx + np.linspace(dx/2, self.box[0] - (dx/2), npores)
            row_y = i*self.unit_cell[1, 1]*dx + (dx/2)*self.unit_cell[1, 1]
            xy_pore_centers[i*npores:(i + 1)*npores, 0] = row_x
            xy_pore_centers[i*npores:(i + 1)*npores, 1] = row_y

        z_separation = self.box[2] / npoints
        column = np.linspace(0, z_separation * (npoints - 1), npoints)

        # For adding noise to quenched disordered configuration
        # columns = np.zeros([npores**2*ncol_per_pore, column.size])
        # xy_noise = np.zeros([npores**2*ncol_per_pore, column.size, 2])
        # shifts = np.zeros([npores**2*ncol_per_pore])
        # shift_range = 0  # fraction of layer to allow random displacement
        # thetas = np.random.uniform(0, 2*np.pi, size=npores**2)  # randomly rotate each pore about the z-axis
        # for i in range(npores**2*ncol_per_pore):
        #     columns[i, :] = z_correlation(column, 20, v=1.2)
        #     xy_noise[i, :, :] = np.random.normal(scale=2.3, size=(column.size, 2))
        #     shifts[i] = shift_range * (z_separation / 2) * np.random.uniform(-1, 1)  # shift column by a random amount

        print('z-spacing: %.2f' % z_separation)
        print('Pore center spacing: %.2f' % dx)

        #thetas = np.random.uniform(0, 2*np.pi, size=npores**2)  # randomly rotate each pore about the z-axis

        for t in range(self.nframes):
            for c in range(npores**2):
                # for each column, choose a random point on the circle with radius, r, centered at the pore center
                # place a column on that point. Equally space remaining columns on circle with reference to that point
                start_theta = np.random.uniform(0, 360) * (np.pi / 180)  # random angle
                # start_theta = thetas[c]
                # start_theta = 0
                theta = 2 * np.pi / ncol_per_pore  # angle between columns
                for a in range(ncol_per_pore):
                    x = r*np.cos(start_theta + a*theta)
                    y = r*np.sin(start_theta + a*theta)
                    start = ncol_per_pore * c * npoints + a * npoints
                    end = ncol_per_pore * c * npoints + (a + 1) * npoints
                    self.locations[t, start:end, :2] = xy_pore_centers[c, :] + [x, y]
                    if noise:
                        shift_range = .66  # fraction of layer to allow random displacements
                        shift = shift_range * (z_separation / 2) * np.random.uniform(-1, 1)  # shift column by a random amount
                    else:
                        shift = 0

                    x_disorder = r*np.random.normal(scale=thermal_disorder[0], size=(end - start))
                    y_disorder = r*np.random.normal(scale=thermal_disorder[1], size=(end - start))
                    z_disorder = z_separation*np.random.normal(scale=thermal_disorder[2], size=(end - start))

                    disorder = np.vstack((x_disorder, y_disorder, z_disorder)).T
                    self.locations[t, start:end, 2] = z_correlation(column, 20, v=1.2) + shift
                    self.locations[t, start:end, :2] += np.random.normal(scale=2.3, size=(column.size, 2))
                    # self.locations[t, start:end, 2] = column + shift

                    # for noise about initial configuration
                    # self.locations[t, start:end, 2] = columns[c*ncol_per_pore + a] + shifts[c*ncol_per_pore + a]
                    # self.locations[t, start:end, :2] += xy_noise[c*ncol_per_pore + a, ...]

                    self.locations[t, start:end, :] += disorder

        from LLC_Membranes.llclib import file_rw

        gamma = 2 * np.pi / 3
        a, b, c = self.box
        A = np.array([a/10, 0, 0])  # vector in x direction
        B = np.array([b/10 * np.cos(gamma), b/10 * np.sin(gamma), 0])  # vector in y direction
        C = np.array([0, 0, c/10])

        unitcell_vectors = np.zeros([self.nframes, 3, 3])
        for i in range(frames):
            # vectors don't change but need them as a trajectory
            unitcell_vectors[i, 0, :] = A
            unitcell_vectors[i, 1, :] = B
            unitcell_vectors[i, 2, :] = C

        file_rw.write_gro_pos(self.locations[-1, ...]/10, 'test.gro', ucell=unitcell_vectors[-1, ...])
        traj = md.formats.TRRTrajectoryFile('test.trr', mode='w', force_overwrite=True)  # create mdtraj TRR trajectory object
        time = np.linspace(0, 1000, self.nframes)  # arbitrary times. Times are required by mdtraj
        traj.write(self.locations/10, time=time, box=unitcell_vectors)  # write the trajectory in .trr format

    def random_layer_rotations(self, npores, ncol_per_pore, r, nlayers, frames=1, thermal_disorder=[0, 0, 0]):

        # create columns
        self.hexagonal_column_grid(npores, ncol_per_pore, r, nlayers, frames=frames, noise=False,
                                   thermal_disorder=thermal_disorder)

        pts_per_pore = nlayers*ncol_per_pore

        for f in range(frames):
            for p in range(npores**2):  # there are npores x npores total pores in the unit cell
                for i in range(nlayers):
                    theta = np.random.uniform(0, 360) * (np.pi / 180)  # random angle by which to rotate layer
                    start = p*(ncol_per_pore*nlayers) + i
                    end = (p + 1) * (ncol_per_pore * nlayers)
                    layer = self.locations[f, start:end:nlayers, :]
                    self.locations[f, start:end:nlayers, :] = rotate_about_z(theta, layer)

    def put_in_box(self):

        zv = [0.0, 0.0, 0.0]  # zero vector
        L = self.box

        # put all atoms inside box - works for single frame and multiframe
        for it in range(self.locations.shape[0]):  # looped to save memory
            self.locations[it, ...] = np.where(self.locations[it, ...] < L, self.locations[it, ...],
                                          self.locations[it, ...] - L)  # get positions in periodic cell
            self.locations[it, ...] = np.where(self.locations[it, ...] > zv, self.locations[it, ...],
                                               self.locations[it, ...] + L)
        #
        # for t in range(self.nframes):
        #     for i, z in enumerate(self.locations[t, :, 2]):
        #         if z < 0:
        #             self.locations[t, i, 2] += self.box[2]
        #         elif z > self.box[2]:
        #             self.locations[t, i, 2] -= self.box[2]

    def compute_structure_factor(self, grid, hexagonal=False, weights=None):

        if hexagonal:
            print("Transforming coordinates to cubic cell")
            self.locations[..., 1] /= np.sin(self.theta)
            self.locations[..., 0] -= self.locations[..., 1] * np.cos(self.theta)

        self.put_in_box()

        # put locations into discrete bins
        # define bin edges in each dimension
        x = np.linspace(0, self.box[0], grid[0] + 1)
        y = np.linspace(0, self.box[1], grid[1] + 1)
        z = np.linspace(0, self.box[2], grid[2] + 1)

        print('Histogramming...')
        H = np.zeros([self.nframes, grid[0], grid[1], grid[2]])
        for f in tqdm.tqdm(range(self.nframes)):
            H[f, ...] = np.histogramdd(self.locations[f, ...], bins=(x, y, z))[0]

        # fourier transform
        print('Computing Fourier Transforms')
        sf = np.zeros([grid[0], grid[1], grid[2]])
        rpi = np.zeros([self.nframes])
        for f in tqdm.tqdm(range(self.nframes)):
            #fft = np.fft.fftn(H[f, ...])
            fft = np.fft.fftn(H[f, ...] - H[f, ...].mean())
            sf += (fft * fft.conjugate()).real

        sf /= (self.nframes * self.locations.shape[1])

        # fft frequencies organized so 0 frequency is at the center in all dimensions
        freq_x = np.fft.fftfreq(grid[0], d=x[1]-x[0])
        ndx = np.argsort(freq_x)
        self.freq_x = freq_x[ndx] * 2 * np.pi

        freq_y = np.fft.fftfreq(grid[1], d=y[1]-y[0])
        ndy = np.argsort(freq_y)
        self.freq_y = freq_y[ndy] * 2 * np.pi

        freq_z = np.fft.fftfreq(grid[2], d=z[1]-z[0])
        ndz = np.argsort(freq_z)
        self.freq_z = freq_z[ndz] * 2 * np.pi

        # reorganize grid
        sf = sf[ndx, :, :]
        sf = sf[:, ndy, :]
        self.sf = sf[:, :, ndz]

        # if hexagonal:
        #
        #     a1 = self.unit_cell[0, :]
        #     a2 = self.unit_cell[1, :]
        #     a3 = self.unit_cell[2, :]
        #
        #     b1 = (np.cross(a2, a3)) / (np.dot(a1, np.cross(a2, a3)))
        #     b2 = (np.cross(a3, a1)) / (np.dot(a2, np.cross(a3, a1)))
        #     b3 = (np.cross(a1, a2)) / (np.dot(a3, np.cross(a1, a2)))
        #
        #     b_inv = np.linalg.inv(np.vstack((b1, b2, b3)))
        #
        #     freq_X, freq_Y, freq_Z = np.meshgrid(self.freq_x, self.freq_y, self.freq_z)
        #     freqs = np.array([freq_X.flatten(), freq_Y.flatten(), freq_Z.flatten(), self.fft.flatten()]).T
        #     freqs[:, :3] = np.matmul(freqs[:, :3], b_inv)
        #     self.fft, edges = np.histogramdd(freqs[:, :3], bins=grid, weights=freqs[:, -1])
        #
        #     self.freq_x = np.array([(edges[0][i - 1] + edges[0][i])/2 for i in range(1, len(edges[0]))])
        #     self.freq_y = np.array([(edges[1][i - 1] + edges[1][i])/2 for i in range(1, len(edges[1]))])
        #     self.freq_z = np.array([(edges[2][i - 1] + edges[2][i])/2 for i in range(1, len(edges[2]))])

    def plot_sf_slice(self, axis, loc, show=False):
        """ Plot 1D slice of structure factor
        :param axis : x, y or z
        :param loc : q values where slice is located perpendicularly. For example, if you take a z-slice then you need
        the qx and qy values that the slice passes through. Pass them in alphabetical order. For example, if you want
        a slice through y, pass in qx then qz. If you want x then qy then qz etc.
        """

        plt.figure()
        axes = {'x': 0, 'y': 1, 'z': 2}

        if axis == 'x':
            ndx_ax1 = np.argmin(np.abs(self.freq_y - loc[0]))
            ndx_ax2 = np.argmin(np.abs(self.freq_z - loc[1]))
            self.slice = self.sf[:, ndx_ax1, ndx_ax2]
            plt.plot(self.freq_x, self.slice)
        elif axis == 'y':
            ndx_ax1 = np.argmin(np.abs(self.freq_x - loc[0]))
            ndx_ax2 = np.argmin(np.abs(self.freq_z - loc[1]))
            self.slice = self.sf[ndx_ax1, :, ndx_ax2]
            plt.plot(self.freq_y, self.slice)
        elif axis == 'z':
            ndx_ax1 = np.argmin(np.abs(self.freq_x - loc[0]))
            ndx_ax2 = np.argmin(np.abs(self.freq_y - loc[1]))
            self.slice = self.sf[ndx_ax1, ndx_ax2, :]
            plt.plot(self.freq_z, self.slice)
        else:
            print('invalid axis chosen for slice of structure factor')

        plt.title('%s slice' % axis)
        np.savez_compressed('correlation.npz', freq_z=self.freq_z, slice=self.slice)

        if show:
            plt.show()

    def scatter3d(self, show=True):
        """
        Create a 3D scatter plot of data
        :param data: x, y, z value to be plotted - np.array([npts, 3])
        :param colorbar: whether to include a colorbar (not implemented)
        :param show: whether to show the plot immediately
        :return: n/a
        """

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(self.locations[0, :, 0], self.locations[0, :, 1], self.locations[0, :, 2])
        plt.xlabel('x')
        plt.ylabel('y')
        ax.set_zlabel('z')

        if show:
            plt.show()

    def angle_average(self, ucell=None, NBR=80, rmax=-1, zbins=-1, zmax=-1, plot=True, show=False, save=False):

        ES = RegularGridInterpolator((self.freq_x, self.freq_y, self.freq_z), self.sf, bounds_error=False)

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
            MAX = np.amax(self.angle_averaged)*.5
            lvls = np.linspace(MIN, MAX, 200)
            ax.contourf(self.r_angle_averaged, self.z_angle_averaged, self.angle_averaged.T, levels=lvls, cmap='jet',
                            extend='max')
            plt.xlabel('$q_r (\AA^{-1})$')
            plt.ylabel('$q_z (\AA^{-1})$')

        if save:
            plt.savefig('rzplot.png')
        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    tol = 0.0001
    box = [float(x) for x in args.box]
    grid = [int(x) for x in args.grid]
    thermal_disorder = [float(x) for x in args.thermal_disorder]

    t = Trajectory()

    bounds = [[0, box[0]], [0, box[1]], [0, box[2]]]

    if 1 - (box[2] % args.dbwl) / args.dbwl > tol:
        dbwl = box[2] / int(box[2] / args.dbwl)
        print('WARNING: the chosen z-spacing will result in a disruption of periodicity. The z-spacing has been '
              'rescaled to %.2f as a fix. Choose a different spacing or different box dimensions if this is not '
              'sufficient' % dbwl)
    else:
        dbwl = args.dbwl

    if args.random_columns:

        t.box = box

        if args.hexagonal:

            t.set_up_hexagonal(args.cell_theta)
            t.hexagonal_column_grid(args.npores, args.ncol_per_pore, args.pore_radius, int(float(args.box[2]) / dbwl),
                                    frames=args.nframes, noise=args.nonoise, thermal_disorder=thermal_disorder)
            t.atomic_form_factor = np.ones(t.locations.shape[1])
        else:

            t.square_column_grid(args.ncolumns, int(float(args.box[2]) / dbwl), frames=args.nframes, bounds=bounds,
                                 noise=args.nonoise)

    elif args.random_layers:

            t.box = box
            t.set_up_hexagonal(args.cell_theta)
            t.random_layer_rotations(args.npores, args.ncol_per_pore, args.pore_radius, int(float(args.box[2]) / dbwl),
                                     frames=args.nframes, thermal_disorder=thermal_disorder)
    else:

        print('There is no structure of which to calculate the structure factor. Please create a custom configuration'
              'using arguments or pass in a GROMACS .gro or trajectory file (.xtc or .trr).')
        exit()

    # plot points in 3D before any modification
    t.scatter3d(show=False)

    t.compute_structure_factor(grid, hexagonal=args.hexagonal)
    t.plot_sf_slice('z', [0, 0], show=False)
    t.angle_average(plot=True, show=False, save=True)

    rpi_index = np.argmin(np.abs(t.freq_z + (2*np.pi/dbwl)))
    print('R-pi intensity: %.2f' % np.amax(t.slice[(rpi_index - 1): (rpi_index + 1)]))

    # fit lorentzian to R-pi
    t.plot_sf_slice('y', [0, 1.7], show=False)

    #np.savez_compressed('perfect_100pores.npz', freq_y=t.freq_y, slice=t.slice)

    # qbound = 0.5  # distance from qz axis to check for peaks
    # lower = np.argmin(np.abs(t.freq_y + qbound))
    # upper = np.argmin(np.abs(t.freq_y - qbound))
    # upper += 1
    #
    # peaks = find_peaks(t.freq_y[lower:upper], t.sf[np.argmin(np.abs(t.freq_x)), lower:upper, rpi_index], tol=5)
    # peaks += lower
    #
    # #peaks = np.linspace(0, t.freq_y.size - 1, t.freq_y.size, dtype=int)  # for disordered columns
    #
    # plt.scatter(t.freq_y[peaks], t.sf[np.argmin(np.abs(t.freq_x)), peaks, rpi_index])
    #
    # # Lorentzian fit (not as good as gaussian)
    # # p = np.array([0.1, 0, t.locations.shape[1]])
    # # solp_lorentz, cov_x = curve_fit(lorentz, t.freq_y[peaks], t.sf[np.argmin(np.abs(t.freq_x)), peaks, rpi_index], p,
    # #                         bounds=[[0, -np.inf, 0], [np.inf, np.inf, np.inf]])
    # #
    # # plt.plot(t.freq_y, lorentz(t.freq_y, solp_lorentz[0], solp_lorentz[1], solp_lorentz[2]), '--', color='red',
    # #          label='Lorentzian', linewidth=2)
    # #
    # # print("Lorentzian FWHM = %.2f A^-1" % solp_lorentz[0])
    #
    # p = np.array([0, 0.3, t.locations.shape[1], 1])
    # solp, cov_x = curve_fit(gaussian, t.freq_y[peaks], t.sf[np.argmin(np.abs(t.freq_x)), peaks, rpi_index], p,
    #                         bounds=([-np.inf, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf]))
    #
    # plt.plot(t.freq_y, gaussian(t.freq_y, solp[0], solp[1], solp[2], solp[3]), '--', color='green', label='Gaussian',
    #          linewidth=2)
    #
    # print("Gaussian FWHM = %.3f +/- %.3f A^-1" % (2*np.sqrt(2*np.log(2))*solp[1],
    #                                        2 * np.sqrt(2 * np.log(2)) * cov_x[1, 1] ** 0.5))
    plt.legend()
    plt.xlabel('$q_y (\AA^{-1}$)')
    plt.ylabel('Intensity')
    plt.tight_layout()

    plt.show()