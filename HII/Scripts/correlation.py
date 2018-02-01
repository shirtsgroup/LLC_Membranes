#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import Atom_props
import tqdm
import matplotlib.pyplot as plt
from matplotlib import ticker
from llclib import physical
from llclib import transform
from llclib import fast_rotate
from place_solutes import trace_pores
from scipy import ndimage


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file. Make sure to '
                                                                            'preprocess with gmx trjconv -pbc whole')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-a', '--atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Name of atoms to calculate'
                        'correlation function with respect to. The center of mass will be used')
    parser.add_argument('--itp', default='/home/bcoscia/PycharmProjects/GitHub/HII/top/Monomer_Tops/NAcarb11V.itp')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Start frame')
    parser.add_argument('-bins', nargs='+', default=100, type=int, help='Integer or array of bin values. If more than'
                        'one value is used, order the inputs according to the order given in args.axis')
    parser.add_argument('-m', '--monomers_per_layer', default=5, type=int, help='Number of monomers per layer')
    parser.add_argument('-s', '--slice', default='xy', help='Slice to be visualized')
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

    args = parser.parse_args()

    return args


def com(pos, mass):
    """
    Calculate center of mass of groups of atoms. Assumes groups are sequentially numbered.
    :param pos: Positions of all atoms for all frames in mdtraj format
    :param mass: mass of each atom in the group whose center of mass will be calculated
    :return: trajectory of center of mass coordinates
    """

    n = len(mass)  # number of atoms
    nT = pos.shape[0]  # number of frames
    ncom = int(pos.shape[1]/n)  # number of centers of mass to calculate at each frame
    mres = np.sum(mass)

    centers = np.zeros([pos.shape[0], ncom, 3])  # will hold positions of all centers of masses

    for f in range(pos.shape[0]):  # loop through all trajectory frames
        for i in range(ncom):
            w = (pos[f, i*n:(i+1)*n, :].T * mass).T  # weight each atom in the residue by its mass
            centers[f, i, :] = np.sum(w, axis=0) / mres  # sum the coordinates and divide by the mass of the residue

    return centers


def duplicate_periodically(pts, box):
    """
    Duplicate points periodically in the +/- xyz directions. Assumes y box vector is angled with respect to x and
    the z vector points straight up
    :param pts: pts to be duplicated
    :param box: unitcell vectors in mdtraj format (use t.unitcell_vectors)
    :return: periodically extended system
    """

    nT = pts.shape[0]
    npts = pts.shape[1]
    p = np.zeros([nT, npts*27, 3])

    p[:, :npts, :] = pts

    for t in range(nT):

        p[t, npts:2*npts, :] = pts[t, :, :] + box[t, 2, :]
        p[t, 2*npts:3*npts, :] = pts[t, :, :] - box[t, 2, :]

        for i in range(3, 6):
            p[t, i*npts:(i+1)*npts, :] = p[t, (i-3)*npts:(i-2)*npts] + box[t, 0, :]

        for i in range(6, 9):
            p[t, i*npts:(i+1)*npts, :] = p[t, (i-6)*npts:(i-5)*npts] - box[t, 0, :]

        for i in range(9, 18):
            p[t, i*npts:(i+1)*npts, :] = p[t, (i-9)*npts:(i-8)*npts] + box[t, 1, :]

        for i in range(18, 27):
            p[t, i*npts:(i+1)*npts, :] = p[t, (i-18)*npts:(i-17)*npts] - box[t, 1, :]

    return p

if __name__ == "__main__":

    args = initialize()

    npores = 4

    dimensions = []
    for i in args.slice:
        if i == 'x':
            dimensions.append(0)
        elif i == 'y':
            dimensions.append(1)
        elif i == 'z':
            dimensions.append(2)

    ndimensions = 3

    if type(args.bins) is list:
        if len(args.bins) > 1:
            bins = np.array(args.bins)
        else:
            bins = np.array([args.bins[0]]*ndimensions)
    else:
        bins = np.array([args.bins]*ndimensions)

    t = md.load(args.traj, top=args.gro)[args.begin:]
    print('Trajectory loaded')

    mass = [Atom_props.mass[i] for i in args.atoms]  # mass of reference atoms

    L = np.zeros([3])  # average box vectors in each dimension
    for i in range(3):
        L[i] = np.mean(np.linalg.norm(t.unitcell_vectors[:, i, :], axis=1))

    if args.range:
        hist_range = []
        for i in range(ndimensions):
            hist_range.append([])
            hist_range[i].append(float(args.range[i*2]))
            hist_range[i].append(float(args.range[i*2 + 1]))
    else:
        hist_range = []
        for i in range(ndimensions):
            hist_range.append([-L[i]/2, L[i]/2])

    keep = [a.index for a in t.topology.atoms if a.name in args.atoms]  # indices of atoms to keep

    pore_spline = np.zeros([t.n_frames, npores*args.layers, 3])
    for frame in range(t.n_frames):
        pore_spline[frame, :, :] += trace_pores(t.xyz[frame, keep, :], t.unitcell_vectors[frame, :2, :2], args.layers)

    ###################### 3D center of mass #########################

    keep = [a.index for a in t.topology.atoms if a.name in args.atoms]  # indices of atoms to keep

    natoms = len(args.atoms)
    monomers_per_layer = int(len(keep) / args.layers / npores / len(args.atoms))  # divide by len(args.atoms) because of com

    center_of_mass = com(t.xyz[:, keep, :], mass)  # calculate centers of mass of atom groups

    com_per_pore = int(center_of_mass.shape[1] / npores)

    periodic_pts = duplicate_periodically(center_of_mass, t.unitcell_vectors)

    correlation = np.zeros(bins)

    x = np.array([1, 0, 0])
    xnorm = np.linalg.norm(x)

    if args.load:

        corr = np.load('correlation_%s%s.npz' %(args.slice[0], args.slice[1]))
        correlation = corr['correlation']
        edges = corr['edges']
        frames = corr['frames']

    else:
        frames = t.n_frames

        if args.slice == 'xy' or args.slice == 'yx':

            for frame in tqdm.tqdm(range(frames)):
                for p in range(npores):
                    for l in range(args.layers):
                        for a in range(monomers_per_layer):

                            pt = p*com_per_pore + l*monomers_per_layer + a  # index of center of mass reference point
                            point = periodic_pts[frame, pt, :]  # coordinates of center of mass reference point
                            v = pore_spline[frame, p*args.layers + l, :] - point  # vector from point to pore center
                            rot = fast_rotate.rotate_vector(periodic_pts[frame, :, :], v, x)  # rotate all points by angle
                            trans = rot - rot[pt, :]  # translate all point so reference point is at the center
                            H, edges = np.histogramdd(trans, bins=bins, range=hist_range)
                            correlation += H
        else:

            for frame in tqdm.tqdm(range(frames)):
                for p in range(npores):
                    for l in range(args.layers):
                        for a in range(monomers_per_layer):

                            pt = p*com_per_pore + l*monomers_per_layer + a  # index of center of mass reference point
                            translated = periodic_pts[frame, :, :] - periodic_pts[frame, pt, :]  # make pt the origin
                            H, edges = np.histogramdd(translated, bins=bins, range=hist_range)
                            correlation += H

        np.savez_compressed('correlation_%s%s' %(args.slice[0], args.slice[1]), correlation=correlation, edges=edges, frames=frames)

    normalization = frames * center_of_mass.shape[1]**2 / np.prod(bins)
    correlation /= normalization

    centers1 = [edges[dimensions[0]][i] + ((edges[dimensions[0]][i + 1] - edges[dimensions[0]][i])/2) for i in range(bins[0])]
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

    exit()

    ################### 1D Center of mass method ###################
    # keep = [a.index for a in t.topology.atoms if a.name in args.atoms]  # indices of atoms to keep
    #
    # center_of_mass = com(t.xyz[:, keep, :], mass)  # calculate centers of mass of atom groups
    #
    # periodic_pts = duplicate_periodically(center_of_mass, t.unitcell_vectors)
    #
    # z = np.zeros([args.bins])
    #
    # for frame in tqdm.tqdm(range(t.n_frames)):
    #     for i in range(400):
    #         H, edges = np.histogram(periodic_pts[frame, i, 2] - periodic_pts[frame, :, 2], bins=args.bins, range=(0, 4))
    #         # H, edges = np.histogramdd(center_of_mass[i, :, :], bins=bins)
    #         z += H
    #
    # centers = [edges[i] + ((edges[i + 1] - edges[i])/2) for i in range(args.bins)]
    #
    # plt.plot(centers, z)
    # plt.show()

    ##################################################################

    ###################### 2D center of mass #########################
    ######## specialized for rings surrounding pore center ###########
    # keep = [a.index for a in t.topology.atoms if a.name in args.atoms]  # indices of atoms to keep
    #
    # natoms = len(args.atoms)
    # monomers_per_layer = int(len(keep) / args.layers / npores / len(args.atoms))  # divide by len(args.atoms) because of com
    #
    # center_of_mass = com(t.xyz[:, keep, :], mass)  # calculate centers of mass of atom groups
    #
    # com_per_pore = int(center_of_mass.shape[1] / npores)
    #
    # periodic_pts = duplicate_periodically(center_of_mass, t.unitcell_vectors)
    #
    # z = np.zeros(bins)
    #
    # x = np.array([1, 0, 0])
    # xnorm = np.linalg.norm(x)
    #
    # if args.axis == 'xy' or args.axis == 'yx':
    #
    #     for frame in tqdm.tqdm(range(t.n_frames)):
    #         for p in range(npores):
    #             for l in range(args.layers):
    #                 for a in range(monomers_per_layer):
    #
    #                     pt = p*com_per_pore + l*monomers_per_layer + a  # index of center of mass reference point
    #                     point = periodic_pts[frame, pt, :]  # coordinates of center of mass reference point
    #                     v = pore_spline[frame, p*args.layers + l, :] - point  # vector from point to pore center
    #                     rot = transform.rotate_vector(periodic_pts[frame, :, :], v, x)  # rotate all points by angle
    #                     trans = rot - rot[pt, :]  # translate all point so reference point is at the center
    #                     H, xedges, yedges = np.histogram2d(trans[:, 0], trans[:, 1], bins=bins, range=hist_range)
    #                     z += H
    # else:
    #
    #     for frame in tqdm.tqdm(range(t.n_frames)):
    #         for p in range(npores):
    #             for l in range(args.layers):
    #                 for a in range(monomers_per_layer):
    #
    #                     pt = p*com_per_pore + l*monomers_per_layer + a  # index of center of mass reference point
    #                     translated = periodic_pts[frame, :, :] - periodic_pts[frame, pt, :]  # make pt the center
    #                     H, xedges, yedges = np.histogram2d(translated[:, dimensions[0]], translated[:, dimensions[1]],
    #                                                        bins=bins, range=hist_range)
    #
    #                     z += H
    #
    # z /= (t.n_frames*npores*args.layers*monomers_per_layer)  # normalize although this does nothing to the visualization
    # # z /= t.n_frames*npores
    #
    # xcenters = [xedges[i] + ((xedges[i + 1] - xedges[i])/2) for i in range(bins[0])]
    # ycenters = [yedges[i] + ((yedges[i + 1] - yedges[i])/2) for i in range(bins[1])]
    #
    # # remove center region of 2D histogram since it is brightest and contains no useful information
    # # first find the center bins (index)
    #
    # for x in range(len(xcenters)):
    #     for y in range(len(ycenters)):
    #         if np.linalg.norm([xcenters[x], ycenters[y]]) < args.block_radius:
    #             z[x, y] = 0
    #
    # print(np.mean(z))
    # # Imax = np.amax(z)
    #
    # # heatmap = plt.imshow(z.T/(args.scale*Imax), extent=[xcenters[0], xcenters[-1], ycenters[0], ycenters[-1]], cmap='jet')
    # # plt.imshow(z.T/Imax, extent=[xcenters[0], xcenters[-1], ycenters[0], ycenters[-1]], cmap='jet')
    #
    # heatmap = plt.imshow(z.T/np.mean(z), extent=[xcenters[0], xcenters[-1], ycenters[0], ycenters[-1]], cmap='jet')
    # plt.imshow(z.T/np.mean(z), extent=[xcenters[0], xcenters[-1], ycenters[0], ycenters[-1]], cmap='jet')
    #
    # cbar = plt.colorbar(heatmap)
    # tick_locator = ticker.MaxNLocator(nbins=5)
    # cbar.locator = tick_locator
    # cbar.update_ticks()
    # # cbar.ax.set_yticklabels(['0', '0.2', '0.4', '0.6', '0.8', '1'])
    # plt.xlabel('%s dimension (nm)' % args.axis[0])
    # plt.ylabel('%s dimension (nm)' % args.axis[1])
    #
    # # plt.plot(centers, z)
    # plt.tight_layout()
    # if args.save:
    #     plt.savefig('%s.png' % args.save)
    # if not args.noshow:
    #     plt.show()
    # exit()
    #####################################################################

    ################### 1D Averaging Method ##########################
    # z = np.zeros([len(args.atoms), 1000])
    #
    # for j, atom in enumerate(args.atoms):
    #
    #     #keep = [a.index for a in t.topology.atoms if a.name in args.atoms]  # indices of atoms to keep
    #     keep = [a.index for a in t.topology.atoms if a.name == atom]  # indices of atoms to keep
    #
    #     # center_of_mass = com(t.xyz[:, keep, :], mass)  # calculate centers of mass of atom groups
    #
    #     periodic_pts = duplicate_periodically(t.xyz[:, keep, :], t.unitcell_vectors)
    #
    #     nbins = args.bins
    #     # bins = [nbins, int(nbins*np.sin(np.pi/3)), 20*nbins]
    #     # g = np.zeros(bins)
    #     # z = np.zeros([len(args.atoms), 1000])
    #
    #     # for t in tqdm.tqdm(range(t.n_frames)):
    #     #     for i in range(400):
    #     #         H, edges = np.histogramdd(periodic_pts[t, :, :] - periodic_pts[t, i, :], bins=bins)
    #     #         # H, edges = np.histogramdd(center_of_mass[i, :, :], bins=bins)
    #     #         g += H
    #     for frame in tqdm.tqdm(range(t.n_frames)):
    #         for i in range(400):
    #             H, edges = np.histogram(periodic_pts[frame, i, 2] - periodic_pts[frame, :400, 2], bins=1000, range=(0, 4))
    #             # H, edges = np.histogramdd(center_of_mass[i, :, :], bins=bins)
    #             z[j, :] += H
    #
    #     z[j, :] /= (400 * t.n_frames)
    #
    # z_avg = np.mean(z, axis=0)
    #
    # centers = [edges[i] + ((edges[i + 1] - edges[i])/2) for i in range(1000)]
    #
    # plt.plot(centers, z_avg)
    # plt.show()
    # exit()
    #####################################################################

    ############################## Animation ############################
    # g /= t.n_frames
    # # fig = plt.figure()
    #
    # from matplotlib import animation
    #
    # fig = plt.figure()
    # hmap = plt.imshow(g[:, 15, :])
    # plt.show()
    #
    # def update(i):
    #
    #     if bins[0] > i:
    #         hmap = plt.imshow(g[i, :, :].T)
    #     elif (bins[1] + bins[0]) > i >= bins[0]:
    #         hmap = plt.imshow(g[:, i - bins[0], :].T)
    #     else:
    #         hmap = plt.imshow(g[:, :, i - bins[0] - bins[1]].T)
    #
    #     return hmap,
    #
    # # call the animator.  blit=True means only re-draw the parts that have changed.
    # anim = animation.FuncAnimation(fig, update, frames=sum(bins), interval=50, blit=True)
    #
    # plt.show()








