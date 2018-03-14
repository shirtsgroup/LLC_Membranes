#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import tqdm
from scipy.interpolate import RegularGridInterpolator


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
                locations[l*nmon_per_layer + i, :] = [r*np.cos(i*theta), r*np.sin(i*theta), l*dbwl]

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

    x_pts = 100
    y_pts = 100
    z_pts = 120
    #grid_pts = 120
    nspikes = 20
    dbwl = 3.7
    xsigma = 0.01
    ysigma = 0.01
    zsigma = 0.01
    offset = True
    rd = False # radial displacement
    nmon_per_layer = 6
    rpore = 5

    x = np.linspace(0, nspikes*dbwl, x_pts + 1)
    y = np.linspace(0, nspikes*dbwl, y_pts + 1)
    z = np.linspace(0, nspikes*dbwl, z_pts + 1)
    X, Y, Z = np.meshgrid(x[:-1], y[:-1], z[:-1])

    xycenter = np.array([x[-10] / 2, y[-12] / 2])

    #locations = pore(nmon_per_layer, rpore, nspikes, dbwl, xycenter, offset=offset, radial_displacement=rd)  # define points exactly
    locations = column(nspikes, dbwl, xycenter)
    # scatter3d(locations)
    # exit()

    H, edges = np.histogramdd(locations, bins=(x, y, z))  # descretize points

    nonzeros = np.nonzero(H)  # find the indices of nonzero values

    pts = np.zeros_like(locations)  # array of points that sit on bin locations
    pts[:, 0] = x[nonzeros[0]]
    pts[:, 1] = y[nonzeros[1]]
    pts[:, 2] = z[nonzeros[2]]

    grid = np.zeros([x_pts, y_pts, z_pts, 4])
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

    # plot 1D cross section through centers of gaussian spheres
    # plt.plot(z[:-1], grid[grid_center, grid_center, :, 3])
    # plt.title('1D z-cross-section of real space points')
    #plt.show()
    #exit()

    # points instead of gaussians
    #grid_center = grid_pts // 2
    #grid[grid_center, grid_center, points, 3] = 1
    #fgrid[grid_center, grid_center, points] = 1

    grid_flat = grid.reshape(x_pts*y_pts*z_pts, 4)
    #nonzero = np.ma.masked_where(grid[:, 3] == 0, grid[:, 3])
    #grid[:, 3] = nonzero # mask zero values for plotting

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
    freq_x = np.fft.fftfreq(fgrid.shape[0], d=x[1]-x[0])
    freq_y = np.fft.fftfreq(fgrid.shape[1], d=y[1]-y[0])
    freq_z = np.fft.fftfreq(fgrid.shape[2], d=zbin)
    ndx = np.argsort(freq_x)
    ndy = np.argsort(freq_y)
    ndz = np.argsort(freq_z)
    freq_x = freq_x[ndx]
    freq_y = freq_y[ndy]
    freq_z = freq_z[ndz]

    fourier_bin = freq_x[1] - freq_x[0]
    print('Fourier spacing : %1.3f' % np.abs(fourier_bin))

    averaged, rfin, zfin = angle_average(freq_x, freq_y, freq_z, fft)

    MIN = np.amin(averaged)
    xlim = .30
    # get max within region that will be plotted
    if rfin[-1] > xlim:
        lower_limit = np.where(freq_x < -xlim)[0][-1] + 1
        upper_limit = np.where(freq_x > xlim)[0][0] - 1
        MAX = np.amax(averaged[lower_limit:upper_limit, lower_limit:upper_limit])
    else:
        MAX = np.amax(averaged)

    NLEVELS = 200
    lvls = np.linspace(MIN, MAX, NLEVELS)  # contour levels

    plt.figure()
    cs = plt.contourf(rfin, zfin[0], averaged.T, levels=lvls, cmap='jet', extend='max')
    plt.xlabel('$q_r$ ($\AA^{-1}$)')
    plt.ylabel('$q_z$ ($\AA^{-1}$)')

    #plt.xlabel('r ' + unitlab, fontsize=14)
    #plt.ylabel('z ' + unitlab, fontsize=14)
    plt.gcf().get_axes()[0].set_ylim(-xlim, xlim)
    plt.gcf().get_axes()[0].set_xlim(-xlim, xlim)
    plt.colorbar(format='%.1f')
    #plt.gcf().get_axes()[0].tick_params(labelsize=14)
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
