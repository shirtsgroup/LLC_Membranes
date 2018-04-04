#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
import tqdm
from scipy.interpolate import RegularGridInterpolator
import argparse
import mdtraj as md


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
    parser.add_argument('-offset_angle', default='symmetric', help='Angle to offset adjacent layer (degrees). Default '
                        'is half of 360/args.nmon')
    parser.add_argument('-rd', '--radial_displaced', help='Distance to radially displaced layers')
    parser.add_argument('-rd_angle', help='Angle (degrees), w.r.t. x-axis, directing where layers will be radially'
                                          'displaced')
    parser.add_argument('-nmon', default=5, type=int, help='Number of monomers per layer (int)')
    parser.add_argument('-pore_radius', type=float, default=5, help='Pore radius (float)')
    parser.add_argument('-box', nargs='+', default=[20, 20, 74], help='length of box vectors. Only orthorhombic boxes'
                                                                      'are implemented')
    parser.add_argument('-stagger', default=False, help='Deviation from uniform spacing of layers (float). e.g. If '
                        'layers are spaced 3.7 apart and you choose -stagger 1, then the first layer would be at zero,'
                        'the second at (3.7 - args.stagger) = 2.7, and the third at 2*3.7=7.4')
    parser.add_argument('-o', '--output', default='rzplot', help='Name of angle averaged plot to save. (no file extension)')
    parser.add_argument('--noshow', action="store_true", help='Do not show plots at end')
    parser.add_argument('-gro', help='Name of coordinate file to fourier transform')
    parser.add_argument('-traj', help='Name of trajectory file to fourier transform')
    parser.add_argument('-avg', action="store_true", help='Valid when -traj is specified. Average histograms of atomic'
                        'positions over trajectory and then take fourier transform once of averaged array')

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


def angle_average(X, Y, Z, SF, ucell=None, NBR=80, rmax=-1, zbins=-1, zmax=-1):

    ES = RegularGridInterpolator((X, Y, Z), SF, bounds_error=False)

    THETA_BINS_PER_INV_ANG = 20.
    MIN_THETA_BINS = 10  # minimum allowed bins
    RBINS = NBR

    if ucell is not None:

        a1 = ucell[0]
        a2 = ucell[1]
        a3 = ucell[2]

        b1 = (np.cross(a2, a3)) / (np.dot(a1, np.cross(a2, a3)))
        b2 = (np.cross(a3, a1)) / (np.dot(a2, np.cross(a3, a1)))
        b3 = (np.cross(a1, a2)) / (np.dot(a3, np.cross(a1, a2)))

        b_inv = np.linalg.inv(np.vstack((b1, b2, b3)))

    if zbins == -1:
        ZBINS = Z.shape[0]  # 400
    else:
        ZBINS = zbins

    if zmax == -1:
        ZMAX = Z[-1]
    else:
        ZMAX = zmax

    XR = (X[-1] - X[0])
    YR = (Y[-1] - Y[0])

    if rmax == -1:
        Rmax = min(XR, YR) / 2.0
        Rmax *= 0.95
    else:
        Rmax = rmax

    rarr, rspace = np.linspace(0.0, Rmax, RBINS, retstep=True)
    #zar = np.linspace(Z[0], Z[-1], ZBINS)
    zar = np.linspace(-ZMAX, ZMAX, ZBINS)

    oa = np.zeros((rarr.shape[0], zar.shape[0]))

    circ = 2.*np.pi*rarr  # circumference

    for ir in range(rarr.shape[0]):

        NTHETABINS = max(int(THETA_BINS_PER_INV_ANG*circ[ir]), MIN_THETA_BINS)  #calculate number of bins at this r
        thetas = np.linspace(0.0, np.pi*2.0, NTHETABINS, endpoint=False)  # generate theta array

        t, r, z = np.meshgrid(thetas, rarr[ir], zar)  # generate grid of cylindrical points

        xar = r*np.cos(t)  # set up x,y coords
        yar = r*np.sin(t)

        pts = np.vstack((xar.ravel(), yar.ravel(), z.ravel())).T  # reshape for interpolation

        if ucell is not None:
            # pts = mc_inv(pts, ucell)
            pts = np.matmul(pts, b_inv)

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


def rescale(coords, dims):
    """
    rescale coordinates so that cell dimensions are constant over the simulation
    returns rescaled coordinates and average length
    """

    avgdims = np.average(dims, axis=0)
    a = avgdims / dims
    rc = coords

    for it in range(rc.shape[0]):
        for i in range(3):
            rc[it, :, i] *= a[it, i]

    # dr = old_div(avgdims,Nspatialgrid)

    return rc, avgdims


def pore(nmon_per_layer, r, layers, dbwl, center, offset=True, offset_angle=0, radial_displacement=False, rd_angle=0, stagger=False):
    """

    :param nmon_per_layer: number of monomers per layer (int)
    :param r: radius of pore (float, angstroms)
    :param layers: number of layers (int)
    :param dbwl : distance between layers (angstroms)
    :param center : xy coordinates of center of pore
    :param offset : True if you want adjacent layers in a parallel displaced configuration. False if you want them sandwiched
    :param stagger: stagger layer spacing by an amount (float)
    :return: locations of all points in space
    """

    locations = np.zeros([nmon_per_layer*layers, 3])

    theta = 2*np.pi / nmon_per_layer  # theta between rings

    displacement = []

    for l in range(layers):
        for i in range(nmon_per_layer):
            if offset or stagger:
                if l % 2 == 0:
                    locations[l*nmon_per_layer + i, :] = [r*np.cos(i*theta), r*np.sin(i*theta), l*dbwl]
                else:
                    if offset:
                        locations[l*nmon_per_layer + i, :] = [r*np.cos(i*theta + offset_angle), r*np.sin(i*theta + offset_angle), l*dbwl]
                    else:
                        locations[l*nmon_per_layer + i, :] = [r*np.cos(i*theta), r*np.sin(i*theta), l*dbwl]
                    if stagger:
                        locations[l*nmon_per_layer + i, 2] = l*dbwl - stagger
            else:
                locations[l*nmon_per_layer + i, :] = [r*np.cos(i*theta + offset_angle), r*np.sin(i*theta + offset_angle), l*dbwl]

            if radial_displacement:
                if l % 2 == 0:
                    displacement.append(l*nmon_per_layer + i)

    locations[:, :2] += center

    if radial_displacement:
        theta = rd_angle
        locations[displacement, :2] += [radial_displacement*np.cos(theta), radial_displacement*np.sin(theta)]

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


def normalize_alkanes(R, Z, Raw_Intensity, inner, outer, angle):
    """
    Plot angular integration of 2D WAXS data bounded by a circle defined by radii 'inner' and 'outer'
    :param R: points in r direction
    :param Z: points in z direction
    :param Raw_Intensity: values at all (R, Z) points on grid
    :param inner: inside radius of region bounding alkane reflections
    :param outer: outside radius of region bounding alkane reflections
    :return: Intensity values normalized by average intensity inside alkane region
    """

    nbins = 45
    bins = np.linspace(-90, 90, nbins)

    bw = 180 / (nbins - 1)

    angles = []
    intensity = []
    for i in range(R.shape[0]):
        for j in range(Z.shape[0]):
            if inner < np.linalg.norm([R[i], Z[j]]) < outer:
                angles.append((180/np.pi)*np.arctan(Z[j]/R[i]))
                intensity.append(Raw_Intensity[i, j])

    inds = np.digitize(angles, bins)

    I = np.zeros([nbins])
    counts = np.zeros([nbins])
    for i in range(len(inds)):
        I[inds[i] - 1] += intensity[i]
        counts[inds[i] - 1] += 1

    #Get average intensity in ring excluding 60 degree slice around top and bottom #######

    bin_range = 180 / nbins  # degrees which a single bin covers

    start = int((angle/2) / bin_range)  # start at the bin which covers -60 degrees and above
    end = nbins - start  # because of symmetry

    total_intensity = np.sum(I[start:end])
    avg_intensity = total_intensity / np.sum(counts[start:end])

    print('Average Intensity in alkane chain region : %s' % avg_intensity)

    return avg_intensity


if __name__ == "__main__":

    # initialize variable based on user inputs
    args = initialize()

    nspikes = args.layers
    dbwl = args.dbwl
    offset = args.offset
    if args.radial_displaced:
        rd = float(args.radial_displaced)
        rd_angle = np.pi * float(args.rd_angle) / 180
    else:
        rd = False
        rd_angle = 0
    nmon_per_layer = args.nmon
    rpore = args.pore_radius
    x_pts = int(args.grid[0])
    stagger = args.stagger
    if stagger:
        stagger = float(stagger)
    offset_angle = args.offset_angle
    if offset_angle == 'symmetric':
        offset_angle = 2*np.pi / nmon_per_layer / 2
    else:
        offset_angle = np.pi * float(offset_angle) / 180

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

        if args.gro:

            # compute structure factor of structural input

            if args.traj:
                t = md.load(args.traj, top=args.gro)
            else:
                t = md.load(args.gro)

            names = [a.element.symbol for a in t.topology.atoms]

            for i in range(len(names)):  # mdtraj interprets sodium wrong
                if names[i] == 'N':
                    names[i] = 'NA'

            import Atom_props
            weights = [Atom_props.radius['%s' % a] for a in names]  # give atoms weights based on atomic radius

            locations = t.xyz

            # convert unit cell -- can be made more general
            theta = 2*np.pi / 3  # hard-coded = bad
            ucell = np.array([[1, 0, 0], [np.cos(theta), np.sin(theta), 0], [0, 0, 1]])

            # convert monoclinic cell to cubic cell
            print("transforming coordinates to monoclinic cell (theta={0:f} deg)".format(theta*180.0/np.pi))
            locations[..., 1] = locations[..., 1] / np.sin(theta)
            locations[..., 0] = locations[..., 0] - locations[..., 1]*np.cos(theta)

            L = np.linalg.norm(t.unitcell_vectors, axis=2)

            locations, L = rescale(locations, L)  # make unit cell constant size -- not sure if this is right to do

            # redefine xyz values of grid
            box = t.unitcell_vectors
            x = np.linspace(0, L[0], int(args.grid[0]))
            y = np.linspace(0, L[1], int(args.grid[1]))  # for converted square box, y box length is same as x
            z = np.linspace(0, L[2], int(args.grid[2]))

            # redefine bins
            xbin = x[1] - x[0]
            ybin = y[1] - y[0]
            zbin = z[1] - z[0]

            zv = [0.0, 0.0, 0.0]  # zero vector

            # put all atoms inside box - works for single frame and multiframe
            for it in range(locations.shape[0]):  # looped to save memory
                locations[it, ...] = np.where(locations[it, ...] < L, locations[it, ...], locations[it, ...] - L)  # get positions in periodic cell
                locations[it, ...] = np.where(locations[it, ...] > zv, locations[it, ...], locations[it, ...] + L)

            # from llclib import file_rw
            # file_rw.write_gro_pos(locations[1, :, :], 'test.gro')
            # exit()
        else:
            # test configurations
            locations = pore(nmon_per_layer, rpore, nspikes, dbwl, xycenter, offset=offset, offset_angle=offset_angle,
                             radial_displacement=rd, rd_angle=rd_angle, stagger=stagger)  # define points exactly

        if args.gro:
            if args.traj:
                if args.avg:
                    H = np.zeros([x.size - 1, y.size - 1, z.size - 1])
                    for f in tqdm.tqdm(range(t.n_frames)):
                        H += np.histogramdd(locations[f, :, :], bins=(x, y, z), weights=weights)[0]
                    H /= t.n_frames
                else:
                    H = np.zeros([t.n_frames, x.size - 1, y.size - 1, z.size - 1])
                    for f in tqdm.tqdm(range(t.n_frames)):
                        H[f, ...] = np.histogramdd(locations[f, :, :], bins=(x, y, z), weights=weights)[0]  # discretize points
            else:
                H, edges = np.histogramdd(locations[0, :, :], bins=(x, y, z), weights=weights)  # discretize points
        else:
            H, edges = np.histogramdd(locations, bins=(x, y, z))

    if not args.gro:

        nonzeros = np.nonzero(H)  # find the indices of nonzero values

        pts = np.zeros_like(locations)  # array of points that sit on bin locations
        if args.dim == 1:
            pts = x[nonzeros[0]]
        if args.dim > 1:
            pts[:, 0] = x[nonzeros[0]]
            pts[:, 1] = y[nonzeros[1]]
            if args.dim > 2:
                pts[:, 2] = z[nonzeros[2]]

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

    if args.gro:

        if args.traj and not args.avg:
            fft = np.zeros([x.size - 1, y.size - 1, z.size - 1])
            for f in tqdm.tqdm(range(t.n_frames)):
                fft += np.abs(np.fft.fftn(H[f, ...]))**2
            fft /= t.n_frames
        else:
            fft = np.abs(np.fft.fftn(H))**2
    else:
        fft = np.abs(np.fft.fftn(fgrid))**2

    freq_x = np.fft.fftfreq(x.size - 1, d=xbin)
    ndx = np.argsort(freq_x)
    freq_x = freq_x[ndx]

    if args.dim > 1:
        freq_y = np.fft.fftfreq(y.size - 1, d=ybin)
        ndy = np.argsort(freq_y)
        freq_y = freq_y[ndy]
        if args.dim > 2:
            freq_z = np.fft.fftfreq(z.size - 1, d=zbin)
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

        # reorganize grid
        fft = fft[ndx, :, :]
        fft = fft[:, ndy, :]
        fft = fft[:, :, ndz]

        if args.gro:
            averaged, rfin, zfin = angle_average(freq_x, freq_y, freq_z, fft, ucell=ucell)
            alkane_intensity = normalize_alkanes(rfin, zfin[0], averaged, 2.23, 2.25, 120)
            MAX = 3.1 * alkane_intensity
            xlim = 4.0
            zlim = 4.0
        else:
            averaged, rfin, zfin = angle_average(freq_x, freq_y, freq_z, fft)
            xlim = 0.4
            zlim = 0.4
            if rfin[-1] > xlim:
                r_lower_limit = np.where(freq_x < -xlim)[0][-1] + 1
                r_upper_limit = np.where(freq_x > xlim)[0][0] - 1
                z_lower_limit = np.where(freq_z < - zlim)[0][-1] + 1
                z_upper_limit = np.where(freq_z > zlim)[0][-1] - 1
                MAX = np.amax(averaged[r_lower_limit:r_upper_limit, z_lower_limit:z_upper_limit])
            else:
                MAX = np.amax(averaged)

        MIN = np.amin(averaged)

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
        plt.savefig('%s.png' % args.output)

        fft_center = fft.shape[0] // 2

        plt.figure()
        plt.plot(freq_z, fft[fft_center, fft_center, :])
        plt.title('z-cross-section of 3D ft')

        plt.figure()
        plt.contourf(freq_x, freq_z, fft[fft_center, :, :].T)
        plt.colorbar()
        plt.title('xz cross-section of 3D ft')

    if not args.noshow:
        plt.show()
