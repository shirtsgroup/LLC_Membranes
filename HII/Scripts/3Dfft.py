#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import tqdm
from scipy.interpolate import RegularGridInterpolator
import argparse


def initialize():

    parser = argparse.ArgumentParser(description='Create configurations in 1, 2 or 3 dimensions and take the discrete'
                                                 'fourier transform of the configuration')

    parser.add_argument('-d', '--dim', default=3, type=int, help='Number of dimensions (int)')
    parser.add_argument('-g', '--grid', nargs='+', default=[100, 100, 100], help='Number of real space grid points in'
                        'each direction (list of ints). Length of array must match number of dimensions')
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of layers in z direction (int)')
    parser.add_argument('-dbwl', default=3.7, type=float, help='Distance between layers (float, angstroms)')
    parser.add_argument('-s', '--sigmas', default=[0.01, 0.01, 0.01], nargs='+', help='Sigma in each direction. (list'
                        'of floats). Length of list must match number of dimensions')
    parser.add_argument('-offset', action="store_true", help='Build offset configuration')
    parser.add_argument('-rd', '--radial_displaced', action="store_true", help='Build radial displaced configuration')
    parser.add_argument('-nmon', default=5, type=int, help='Number of monomers per layer (int)')
    parser.add_argument('-pore_radius', type=float, default=5, help='Pore radius (float)')
    parser.add_argument('-box', nargs='+', default=[20, 20, 74], help='length of box vectors. Only orthorhombic boxes'
                                                                      'are implemented')

    args = parser.parse_args()

    return args


def gaussian1D(x, sigma_x, mu_x):
    c = sigma_x*np.sqrt(2*np.pi)
    return (1/c)*np.exp(-(x-mu_x)**2/(2*sigma_x**2))


def gaussian2D(x, y, sigma_x, sigma_y, mu_x, mu_y):
    c = sigma_x*sigma_y*2*np.pi
    return (1/c)*np.exp(-(x-mu_x)**2/(2*sigma_x**2))*np.exp(-(y-mu_y)**2/(2*sigma_y**2))*np.exp(-(z-mu_z)**2/(2*sigma_z**2))


def gaussian3D(x, y, z, sigma_x, sigma_y, sigma_z, mu_x, mu_y, mu_z):
    c = sigma_x*sigma_y*sigma_z*np.sqrt(2*np.pi)**3

    return (1/c)*np.exp(-(x-mu_x)**2/(2*sigma_x**2))*np.exp(-(y-mu_y)**2/(2*sigma_y**2))*np.exp(-(z-mu_z)**2/(2*sigma_z**2))


def angle_average(X, Y, Z, SF):

    ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

    THETA_BINS_PER_INV_ANG = 20.
    MIN_THETA_BINS = 10  # minimum allowed bins
    RBINS = 400

    ZBINS = Z.shape[0]  # 400

    XR = (X[-1] - X[0])
    YR = (Y[-1] - Y[0])

    Rmax = min(XR, YR) / 2.0
    Rmax *= 0.95

    rarr, rspace = np.linspace(0.0, Rmax, RBINS, retstep=True)
    zar = np.linspace(Z[0], Z[-1], ZBINS)

    oa = np.zeros((rarr.shape[0], zar.shape[0]))
    circ = 2.*np.pi*rarr  # circumference

    for ir in tqdm.tqdm(range(rarr.shape[0])):

        NTHETABINS = max(int(THETA_BINS_PER_INV_ANG*circ[ir]), MIN_THETA_BINS)  #calculate number of bins at this r
        thetas = np.linspace(0.0, np.pi*2.0, NTHETABINS, endpoint=False)  # generate theta array

        t, r, z = np.meshgrid(thetas, rarr[ir], zar)  # generate grid of cylindrical points

        xar = r*np.cos(t)  # set up x,y coords
        yar = r*np.sin(t)

        pts = np.vstack((xar.ravel(), yar.ravel(), z.ravel())).T  # reshape for interpolation

        oa[ir, :] = np.average(ES(pts).reshape(r.shape), axis=1)  # store average values in final array

    mn = np.nanmin(oa)
    oa = np.where(np.isnan(oa), mn, oa)

    rad_avg = np.average(oa)  # ???
    oa /= rad_avg  # normalize

    # set up data for contourf plot by making it symmetrical
    final = np.append(oa[::-1, :], oa[1:], axis=0)  # SF
    rfin = np.append(-rarr[::-1], rarr[1:])  # R
    zfin = np.append(z[:, 0, :], z[1:, 0, :], axis=0)  # Z

    return final, rfin, zfin


def pore(nmon_per_layer, r, layers, dbwl, center, offset=True, radial_displacement=False):
    """

    :param nmon_per_layer: number of monomers per layer (int)
    :param r: radius of pore (float, angstroms)
    :param layers: number of layers (int)
    :param dbwl : distance between layers (angstroms)
    :param center : xy coordinates of center of pore
    :param offset : True if you want adjacent layers in a parallel displaced configuration. False if you want them sandwiched
    :return: locations of all points in space
    """

    locations = np.zeros([nmon_per_layer*layers, 3])

    theta = 2*np.pi / nmon_per_layer  # theta between rings
    half_theta = theta / 2

    displacement = []

    for l in range(layers):
        for i in range(nmon_per_layer):
            if offset:
                if l % 2 == 0:
                    locations[l*nmon_per_layer + i, :] = [r*np.cos(i*theta), r*np.sin(i*theta), l*dbwl]
                else:
                    locations[l*nmon_per_layer + i, :] = [r*np.cos(i*theta + half_theta), r*np.sin(i*theta + half_theta), l*dbwl]
            else:
                locations[l*nmon_per_layer + i, :] = [r*np.cos(i*theta + half_theta), r*np.sin(i*theta + half_theta), l*dbwl]

            if radial_displacement:
                if l % 2 == 0:
                    displacement.append(l*nmon_per_layer + i)

    locations[:, :2] += center

    if radial_displacement:
        r = 11
        theta = 0
        locations[displacement, :2] += [r*np.cos(theta), r*np.sin(theta)]

    return locations


def column(layers, dbwl, center):
    """
    Create a single column of points
    :param layers: number of layers (int)
    :param dbwl: distance between layers (angstroms, float)
    :param center: xy position where column should be located
    :return: points np.array([layers, 3])
    """

    locations = np.zeros([layers, 3])
    z = np.linspace(0, (layers - 1)*dbwl, layers)

    locations[:, :2] = center
    locations[:, 2] = z

    return locations


def scatter3d(data, colorbar=False, show=True):
    """
    Create a 3D scatter plot of data
    :param data: x, y, z value to be plotted - np.array([npts, 3])
    :param colorbar: whether to include a colorbar (not implemented)
    :param show: whether to show the plot immediately
    :return: n/a
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(data[:, 0], data[:, 1], data[:, 2])
    plt.xlabel('x')
    plt.xlabel('y')
    plt.xlabel('z')

    if show:
        plt.show()

if __name__ == "__main__":

    # initialize variable based on user inputs
    args = initialize()

    nspikes = args.layers
    dbwl = args.dbwl
    offset = args.offset
    rd = args.radial_displaced
    nmon_per_layer = args.nmon
    rpore = args.pore_radius
    x_pts = int(args.grid[0])
    if args.dim == 1:
        x = np.linspace(0, nspikes*dbwl, x_pts + 1)
    else:
        x = np.linspace(0, float(args.box[0]), x_pts + 1)
    xbin = x[1] - x[0]
    xsigma = float(args.sigmas[0])
    if args.dim > 1:
        y_pts = int(args.grid[1])
        ysigma = args.sigmas[1]
        # y = np.linspace(0, nspikes*dbwl, y_pts + 1)
        y = np.linspace(0, float(args.box[1]), y_pts + 1)
        ybin = y[1] - y[0]
        xycenter = np.array([x[-1] / 2, y[-1] / 2])
        if args.dim == 2:
            X, Y = np.meshgrid(x[:-1], y[:-1])
        if args.dim > 2:
            z_pts = int(args.grid[2])
            zsigma = args.sigmas[2]
            z = np.linspace(0, nspikes*dbwl, z_pts + 1)
            zbin = z[1] - z[0]
            X, Y, Z = np.meshgrid(x[:-1], y[:-1], z[:-1])

    if args.dim == 1:
        locations = np.linspace(0, (nspikes - 1)*dbwl, nspikes)
        H, edges = np.histogram(locations, bins=x)

    if args.dim == 3:

        locations = pore(nmon_per_layer, rpore, nspikes, dbwl, xycenter, offset=offset, radial_displacement=rd)  # define points exactly
        H, edges = np.histogramdd(locations, bins=(x, y, z))  # discretize points

    # locations = column(nspikes, dbwl, xycenter)
    # scatter3d(locations)
    # exit()

    count = 0
    for i in range(len(H)):
        if H[i] == 0:
            count += 1

    nonzeros = np.nonzero(H)  # find the indices of nonzero values

    pts = np.zeros_like(locations)  # array of points that sit on bin locations
    if args.dim == 1:
        pts = x[nonzeros[0]]
    if args.dim > 1:
        pts[:, 0] = x[nonzeros[0]]
        pts[:, 1] = y[nonzeros[1]]
        if args.dim > 2:
            pts[:, 2] = z[nonzeros[2]]

    print(pts)
    # exit()

    if args.dim == 1:

        fgrid = np.zeros([x_pts])
        for i in tqdm.tqdm(range(pts.shape[0]), unit='points'):
            gauss = gaussian1D(x[:-1], xsigma, locations[i])
            fgrid += gauss

    if args.dim == 2:

        fgrid = np.zeros([x_pts, y_pts])
        for i in tqdm.tqdm(range(pts.shape[0]), unit='points'):
            gauss = gaussian2D(X, Y, xsigma, ysigma, pts[i, 0], pts[i, 1])
            fgrid += gauss

    if args.dim == 3:

        grid = np.zeros([x_pts, y_pts, z_pts, 4])  # for scatter plot
        fgrid = np.zeros([x_pts, y_pts, z_pts])

        dist_bw_pts = z_pts // nspikes
        zbin = dbwl / dist_bw_pts

        for i in range(x_pts):
            for j in range(y_pts):
                for k in range(z_pts):
                    grid[i, j, k, :3] = [x[i], y[j], z[k]]

        for i in tqdm.tqdm(range(pts.shape[0]), unit='points'):
            gauss = gaussian3D(X, Y, Z, xsigma, ysigma, zsigma, pts[i, 0], pts[i, 1], pts[i, 2])
            grid[..., 3] += gauss
            fgrid += gauss

    if args.dim == 1:
        plt.plot(x[:-1], fgrid)

    if args.dim == 2:
        plt.contourf(x[:-1], y[:-1], fgrid)

    if args.dim == 3:

        grid_flat = grid.reshape(x_pts*y_pts*z_pts, 4)
        grid_flat = grid_flat[np.where(grid_flat[:, 3] > 0.01), :]
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        cmap = cm.ScalarMappable(cmap=cm.hsv)
        I = grid_flat[0, :, 3]
        cmap.set_array(I)

        ax.scatter(grid_flat[0, :, 0], grid_flat[0, :, 1], grid_flat[0, :, 2], c=cm.hsv(I/max(I)))
        plt.title('Real space points')
        xy_window = 3*rpore
        ax.set_xlim3d(xycenter[0] - xy_window, xycenter[0] + xy_window)
        ax.set_ylim3d(xycenter[1] - xy_window, xycenter[1] + xy_window)
        fig.colorbar(cmap)

    fft = np.abs(np.fft.fftn(fgrid))**2
    #fft = np.fft.rfftn(fgrid)

    freq_x = np.fft.fftfreq(fgrid.shape[0], d=xbin)
    ndx = np.argsort(freq_x)
    freq_x = freq_x[ndx]

    if args.dim > 1:
        freq_y = np.fft.fftfreq(fgrid.shape[1], d=ybin)
        ndy = np.argsort(freq_y)
        freq_y = freq_y[ndy]
        if args.dim > 2:
            freq_z = np.fft.fftfreq(fgrid.shape[2], d=zbin)
            ndz = np.argsort(freq_z)
            freq_z = freq_z[ndz]

    fourier_bin = freq_x[1] - freq_x[0]
    print('Fourier spacing : %s' % np.abs(fourier_bin))

    if args.dim == 1:

        plt.figure()
        #plt.plot(freq_x[-fft.shape[0]:], fft)
        plt.plot(freq_x, fft)

    if args.dim == 2:

        plt.contourf(freq_x, freq_y, fgrid)

    if args.dim == 3:

        averaged, rfin, zfin = angle_average(freq_x, freq_y, freq_z, fft)

        MIN = np.amin(averaged)
        xlim = .30
        zlim = .30
        # get max within region that will be plotted
        if rfin[-1] > xlim:
            r_lower_limit = np.where(freq_x < -xlim)[0][-1] + 1
            r_upper_limit = np.where(freq_x > xlim)[0][0] - 1
            z_lower_limit = np.where(freq_z < - zlim)[0][-1] + 1
            z_upper_limit = np.where(freq_z > zlim)[0][-1] - 1
            MAX = np.amax(averaged[r_lower_limit:r_upper_limit, z_lower_limit:z_upper_limit])
        else:
            MAX = np.amax(averaged)

        NLEVELS = 200
        lvls = np.linspace(MIN, MAX, NLEVELS)  # contour levels

        plt.figure()
        cs = plt.contourf(rfin, zfin[0], averaged.T, levels=lvls, cmap='jet', extend='max')
        plt.xlabel('$q_r$ ($\AA^{-1}$)')
        plt.ylabel('$q_z$ ($\AA^{-1}$)')
        plt.gcf().get_axes()[0].set_ylim(-zlim, zlim)
        plt.gcf().get_axes()[0].set_xlim(-xlim, xlim)
        plt.colorbar(format='%.1f')
        plt.tight_layout()

        fft_center = fft.shape[0] // 2

        plt.figure()
        plt.plot(freq_z, fft[fft_center, fft_center, :])
        plt.title('z-cross-section of 3D ft')

        plt.figure()
        plt.contourf(freq_x, freq_z, fft[fft_center, :, :].T)
        plt.colorbar()
        plt.title('xz cross-section of 3D ft')

    plt.show()
