#! /usr/bin/env python

""" Calculate structure factor of single 3D configurations or trajectories of configurations """

import tqdm
import argparse
import mdtraj as md
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np


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
    parser.add_argument('-box', nargs='+', default=[20, 20, 37], help='length of box vectors. Only orthorhombic boxes'
                                                                      'are implemented')
    parser.add_argument('-nframes', default=100, type=int, help='Number of frames to create')
    parser.add_argument('--random_columns', action="store_true", help='Create a trajectory with columns made of'
                        'equispaced points but are randomly displaced in the z-diretion with respect to other columns')
    parser.add_argument('-ncol', '--ncolumns', default=10, type=int, help='The number of columns in the x and y '
                        'dimensions. The total number of columns in the unit cell will be ncol^2. In the case of '
                        'hexagonal columns, this defines the number of pore centers.')
    parser.add_argument('-dbwl', default=3.7, type=float, help='Distance between layers (float, angstroms)')
    parser.add_argument('-nonoise', action="store_false", help='Turn off random column displacement')
    # hexagonally packed columns (combined with --random_columns flag)
    parser.add_argument('--hexagonal', action="store_true", help='Create a hexagonal packed columns that mimic the HII'
                                                                 'phase')
    parser.add_argument('--ncol_per_pore', default=5, type=int, help='Number of columns surrounding each pore center')
    parser.add_argument('--pore_radius', default=0.5, type=float, help='Distance to place each column from pore center')
    parser.add_argument('--cell_theta', default=120, type=float, help='Angle between vectors defining xy plane of'
                                                                      'monoclinic box')
    parser.add_argument('-npores', default=4, type=int, help='Number of pores in each dimension, similar to ncol')

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


class Trajectory(object):

    def __init__(self):

        # initialize variables
        self.locations = None
        self.box = None
        self.nframes = 0
        self.fft = None
        self.freq_x = None
        self.freq_y = None
        self.freq_z = None
        self.slice = None
        self.unit_cell = None
        self.theta = 0

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

            # make sure everything is within 'box'
            # if bounds:
            #     for i, z in enumerate(self.locations[t, :, 2]):
            #         if z < bounds[2][0]:
            #             self.locations[t, i, 2] += bounds[2][1]
            #         elif z > bounds[2][1]:
            #             self.locations[t, i, 2] -= bounds[2][1]
        self.put_in_box()

    def hexagonal_column_grid(self, npores, ncol_per_pore, r, npoints, frames=1, noise=True):

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

        print('z-spacing: %.2f' % z_separation)
        print('Pore center spacing: %.2f' % dx)

        for t in range(self.nframes):
            for c in range(npores**2):
                # for each column, choose a random point on the circle with radius, r, centered at the pore center
                # place a column on that point. Equally space remaining columns on circle with reference to that point
                start_theta = np.random.uniform(0, 360) * (np.pi / 180)  # random angle
                theta = 2 * np.pi / ncol_per_pore  # angle between columns
                for a in range(ncol_per_pore):
                    x = r*np.cos(start_theta + a*theta)
                    y = r*np.sin(start_theta + a*theta)
                    start = ncol_per_pore * c * npoints + a * npoints
                    end = ncol_per_pore * c * npoints + (a + 1) * npoints
                    self.locations[t, start:end, :2] = xy_pore_centers[c, :] + [x, y]
                    if noise:
                        shift = (z_separation / 2) * np.random.uniform(-1, 1)  # shift column by a random amount
                    else:
                        shift = 0
                    self.locations[t, start:end, 2] = column + shift

        self.put_in_box()

    def put_in_box(self):

        for t in range(self.nframes):
            for i, z in enumerate(self.locations[t, :, 2]):
                if z < 0:
                    self.locations[t, i, 2] += self.box[2]
                elif z > self.box[2]:
                    self.locations[t, i, 2] -= self.box[2]

    def compute_structure_factor(self, grid, hexagonal=False, weights=None):

        if hexagonal:
            print("Transforming coordinates to cubic cell")
            self.locations[..., 1] /= np.sin(self.theta)
            self.locations[..., 0] -= self.locations[..., 1] * np.cos(self.theta)

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
        fft = np.zeros([grid[0], grid[1], grid[2]])
        for f in tqdm.tqdm(range(self.nframes)):
            fft += np.abs(np.fft.fftn(H[f, ...] - H[f, ...].mean()))**2

        fft /= self.nframes

        # reorganize FT so 0 frequency is at the center in all dimensions
        freq_x = np.fft.fftfreq(grid[0], d=x[1]-x[0])
        ndx = np.argsort(freq_x)
        self.freq_x = freq_x[ndx]

        freq_y = np.fft.fftfreq(grid[1], d=y[1]-y[0])
        ndy = np.argsort(freq_y)
        self.freq_y = freq_y[ndy]

        freq_z = np.fft.fftfreq(grid[2], d=z[1]-z[0])
        ndz = np.argsort(freq_z)
        self.freq_z = freq_z[ndz]

        # reorganize grid
        fft = fft[ndx, :, :]
        fft = fft[:, ndy, :]
        self.fft = fft[:, :, ndz]

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
            self.slice = self.fft[:, ndx_ax1, ndx_ax2]
            plt.plot(self.freq_x, self.slice)
        elif axis == 'y':
            ndx_ax1 = np.argmin(np.abs(self.freq_x - loc[0]))
            ndx_ax2 = np.argmin(np.abs(self.freq_z - loc[1]))
            self.slice = self.fft[ndx_ax1, :, ndx_ax2]
            plt.plot(self.freq_y, self.slice)
        elif axis == 'z':
            ndx_ax1 = np.argmin(np.abs(self.freq_x - loc[0]))
            ndx_ax2 = np.argmin(np.abs(self.freq_y - loc[1]))
            self.slice = self.fft[ndx_ax1, ndx_ax2, :]
            plt.plot(self.freq_z, self.slice)
        else:
            print('invalid axis chosen for slice of structure factor')

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


if __name__ == "__main__":

    args = initialize().parse_args()

    tol = 0.0001
    box = [float(x) for x in args.box]
    grid = [int(x) for x in args.grid]

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

            if t.box[0] != t.box[1]:
                print('WARNING: The x and y box lengths must be equal for this script to properly implement hexagonal '
                      'periodicity. Setting y box length equal to x box length.')
                t.box[1] = t.box[0]

            t.theta = args.cell_theta * np.pi / 180.0  # theta for monoclinic unit cell
            t.unit_cell = np.array([[1, 0, 0], [np.cos(t.theta), np.sin(t.theta), 0], [0, 0, 1]])
            t.hexagonal_column_grid(args.npores, args.ncol_per_pore, args.pore_radius, int(float(args.box[2]) / dbwl),
                                    frames=args.nframes, noise=args.nonoise)
        else:

            t.square_column_grid(args.ncolumns, int(float(args.box[2]) / dbwl), frames=args.nframes, bounds=bounds,
                                 noise=args.nonoise)

    else:

        print('There is no structure of which to calculate the structure factor. Please create a custom configuration'
              'using arguments or pass in a GROMACS .gro or trajectory file (.xtc or .trr).')
        exit()

    # plot points in 3D before any modification
    t.scatter3d(show=False)

    t.compute_structure_factor(grid, hexagonal=args.hexagonal)
    t.plot_sf_slice('z', [0, 0], show=False)

    rpi_index = np.argmin(np.abs(t.freq_z + (1/dbwl)))
    print('R-pi intensity: %.2f' % np.amax(t.slice[(rpi_index - 1): (rpi_index + 1)]))

    plt.show()