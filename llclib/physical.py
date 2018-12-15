#! /usr/bin/env python

from __future__ import division
from builtins import range
from LLC_Membranes.llclib import file_rw, transform
import mdtraj as md
import numpy as np
import matplotlib.path as mplPath
import mdtraj as md
from random import randint
from pymbar import timeseries
import tqdm


class region:
    """
    Define a region as an extrusion of a polygon in the z direction
    """

    def __init__(self, sides):
        """
        :param sides: number of sides making up the region in the xy plane
        :return: region
        """
        self.sides = sides

    def xyregion(self, corners):
        """
        :param corners: points defining the corners of the polygon making up the xy region
        :return: a region defined by corners
        """
        path = mplPath.Path(corners)


def thickness(filename, ref_atoms, grid, *traj, **kwargs):
    """
    :param filename: name of .gro file
    :param ref_atoms: atoms which thickness will be based on
    :param traj: trajectory of positions
    :return: trajectory of thicknesses or single thickness based on max/min z coordinate of reference atoms
    """

    if traj:

        traj = np.asarray(traj)[0]  # optional arguments of the form *args need to be convert back to numpy arrays
        nT = traj.shape[0]  # number of trajectory points

        thick = np.zeros([nT])
        z_max = np.zeros([nT])
        z_min = np.zeros([nT])
        thick_std = np.zeros([nT])
        for t in range(nT):
            z_max_t = max(traj[t, :, 2])
            z_min_t = min(traj[t, :, 2])
            thick[t] = z_max_t - z_min_t
            z_max[t] = z_max_t
            z_min[t] = z_min_t

    else:
        f = open(filename, "r")  # .gro file whose positions of Na ions will be read

        a = []  # list to hold lines of file
        for line in f:
            a.append(line)

        line = 0
        while a[line].count('HII') == 0:
            line += 1

        if grid:

            t = md.load(filename)
            pos = t.xyz[0, :, :]  # positions of all atoms

            if kwargs['exclude']:
                keep = [a.index for a in t.topology.atoms if a.residue.name != kwargs['exclude']]
                pos = t.atom_slice(keep).xyz[0, :, :]

            # define boundaries of each grid area
            grid_res = kwargs['grid_res']

            nregions = (grid_res - 1) ** 2
            g = np.zeros([2, grid_res, grid_res])
            dims = a[-1].split()
            xbox = np.linalg.norm([dims[0], dims[3], dims[4]])
            ybox = np.linalg.norm([dims[1], dims[5], dims[6]])
            yangle = np.arctan(float(dims[1])/abs(float(dims[5])))

            for i in range(grid_res):
                g[0, i, :] = np.linspace(0, xbox, grid_res) + (float(i) / grid_res)*float(dims[5])
                g[1, :, i] = np.linspace(0, ybox, grid_res)*np.sin(yangle)

            corners = np.zeros([nregions, 4, 2])
            zmaxes = np.zeros([nregions])
            zmins = np.zeros([nregions])
            thicks = np.zeros([nregions])

            for i in range(grid_res - 1):
                for j in range(grid_res - 1):
                    # define corners of grid region
                    r = i*(grid_res - 1) + j
                    corners[r, 0, :] = g[:, i, j]
                    corners[r, 1, :] = g[:, i + 1, j]
                    corners[r, 2, :] = g[:, i + 1, j + 1]
                    corners[r, 3, :] = g[:, i, j + 1]

                    # create a region using the corners (corners need to be traced in order)
                    path = mplPath.Path(corners[r, :, :])
                    contained = path.contains_points(pos[:, :2])  # check whether each point is in the region
                    z = pos[np.where(contained), 2]  # get the z position of all atoms contained in the region
                    zmaxes[r] = np.max(z)
                    zmins[r] = np.min(z)
                    thicks[r] = zmaxes[r] - zmins[r]

            # bootstrap to get statistics
            nboot = 2000
            vmax = np.zeros([nboot])
            vmin = np.zeros([nboot])

            for i in range(nboot):
                imax = randint(0, nregions - 1)
                imin = randint(0, nregions - 1)
                vmax[i] = zmaxes[imax]
                vmin[i] = zmins[imin]

            z_max = np.mean(vmax)
            z_min = np.mean(vmin)
            thick = np.mean(vmax - vmin)
            thick_std = np.std(vmax - vmin)

        else:
            z = []  # list to hold z positions of all atoms

            while a[line].count('HII') != 0:
                if str.strip(a[line][11:15]) in ref_atoms:
                    z.append(float(a[line][36:44]))
                line += 1

            z_max = max(z)
            z_min = min(z)
            thick = z_max - z_min
            thick_std = 0

    return thick, z_max, z_min, thick_std


def conc(t, comp, b):
    """
    Calculate the concentration of the specified component
    :param t: mdtraj trajectory object for system being studied
    :param comp: component which you want the concentration of
    :param b: buffer. distance into membrane to go before starting calculation
    :return: concentration
    """

    box = t.unitcell_vectors
    equil = timeseries.detectEquilibration(box[:, 2, 2])[0]
    thick = np.mean(box[equil:, 2, 2])

    z_max = thick
    z_min = 0
    buffer = thick*b
    z_max -= buffer
    z_min += buffer
    thick = z_max - z_min

    # Calculate concentration (an average of all frames)
    keep = [a.index for a in t.topology.atoms if a.name == comp]
    t_comp = t.atom_slice(keep)

    pos = t_comp.xyz
    ncomp = pos.shape[1]  # number of components in the simulation which you want the concentration of
    nT = pos.shape[0]

    if b > 0:
        count = np.zeros([nT])
        box_vol = np.zeros([nT])
        cross = np.zeros([nT])
        for t in range(nT):
            x_dim = np.linalg.norm(box[t, 0, :])
            y_dim = np.linalg.norm(box[t, 1, :])
            cross[t] = x_dim*y_dim
            box_vol[t] = x_dim*y_dim*thick
            for c in range(ncomp):
                if z_max >= pos[t, c, 2] >= z_min:
                    count[t] += 1
    else:
        count = ncomp*np.ones([nT])
        box_vol = np.zeros([nT])
        cross = np.zeros([nT])
        for t in range(nT):
            x_dim = np.linalg.norm(box[t, 0, :])
            y_dim = np.linalg.norm(box[t, 1, :])
            cross[t] = x_dim*y_dim
            box_vol[t] = x_dim*y_dim*thick

    factor = 1 / (1*10**-27)  # convert from ions/nm^3 to ions/m^3
    conc = np.zeros([nT])
    for c in range(nT):
        conc[c] = (count[c] / box_vol[c]) * factor

    avg_conc = np.mean(conc)
    std = np.std(conc)
    avg_cross = np.mean(cross)

    return avg_conc, std, avg_cross, thick, z_max, z_min


def avg_pore_loc(npores, pos, box, buffer=0, spline=False, npts=20, progress=False, bins=False):
    """ Calculate average pore location for each pore at each frame

    :param no_pores: the number of pores in the unit cell
    :param pos: the coordinates of the component(s) which you are using to locate the pore centers
    :param box: box vectors (t.unitcell_vectors when trajectory is load with mdtraj)
    :param buffer: fraction (of membrane thickness) of top and bottom of membrane to exclude from p2p calculations
    :param spline: trace pore centers with a spline
    :param npts: number of points making up the spline in each pore
    :param progress: show progress bar while constructing splines
    :param bins: return the bin centers of each spline for plotting purposes

    :type no_pores: int
    :type pos: numpy.ndarray, shape(ncomponents, 3) or numpy.ndarray, shape(nframes, ncomponents, 3)
    :type buffer: float
    :type spline: bool
    :type box: numpy.ndarray, shape(nframes, 3, 3)
    :type npts: int
    :type progress: bool
    :type bins: bool

    :return: numpy array containing the x, y coordinates of the center of each pore at each frame
    """

    # Find the average location of the pores w.r.t. x and y

    if spline:
        if box is None:
            print('You must supply box vectors if you are to trace the pores with a spline')
            exit()
        else:

            print('Calculating pore spline...')
            centers, bin_centers = trace_pores(pos, box, npts, npores=4, progress=progress)

            if bins:
                return centers, bin_centers
            else:
                return centers
    else:

        if len(pos.shape) == 3:  # multiple frames

            nT = np.shape(pos)[0]
            comp_ppore = np.shape(pos)[1] // npores

            p_center = np.zeros([nT, npores, 2])

            for i in range(nT):

                positions = wrap_box(pos[i, ...], box[i, ...])

                if buffer > 0:

                    include = np.full(pos.shape[1], True)

                    include[np.where(pos[i, :, 2] > box[i, 2, 2] + buffer)] = False
                    include[np.where(pos[i, :, 2] < buffer)] = False

                    for j in range(npores):
                        p_center[i, j, :] = positions[comp_ppore * j:comp_ppore * (j + 1), :2].mean(axis=0)
                        count = 0
                        for k in range(comp_ppore * j, comp_ppore * (j + 1)):
                            if include[k]:
                                p_center[i, j, :] += positions[k, :2]
                                count += 1
                        p_center[i, j, :] /= count  # take the average

                else:

                    for j in range(npores):
                        p_center[i, j, :] = positions[comp_ppore*j:comp_ppore*(j + 1), :2].mean(axis=0)

        elif len(pos.shape) == 2:  # single frame

            comp_ppore = pos.shape[0] // npores

            p_center = np.zeros([npores, 2])

            for j in range(npores):
                for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                    p_center[j, :] += pos[k, :2]
                p_center[j, :] /= comp_ppore

        else:
            return 'Please use a position array with valid dimensions'
            exit()

        return p_center


def p2p(p_centers, distances):
    """
    :param p_centers: the x, y locations of the pore centers in the format return from avg_pore_loc()
    :param distances: the number of distinct distances between pores
    :return: all of the pore to pore distances
    """
    nT = np.shape(p_centers)[2]
    p2ps = np.zeros([distances, nT])  # distances in the order 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
    for i in range(nT):
        # So ugly ... sadness :(
        p2ps[0, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 1, i])
        p2ps[1, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 2, i])
        p2ps[2, i] = np.linalg.norm(p_centers[:, 0, i] - p_centers[:, 3, i])
        p2ps[3, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 2, i])
        p2ps[4, i] = np.linalg.norm(p_centers[:, 1, i] - p_centers[:, 3, i])
        p2ps[5, i] = np.linalg.norm(p_centers[:, 2, i] - p_centers[:, 3, i])

    return p2ps


def limits(pos, pcenters):
    """
    Estimate the pore 'radius' based on the position of some component and it's maximum deviation from the pore center
    :param: pos: the positions of all atoms included in making the estimate
    :param: pcenters: the x,y positions of the pore centers for each frame
    :return: an approximate pore radius. Beyond which, we have entered the alkane region
    """

    nT = pcenters.shape[0]
    npores = pcenters.shape[1]
    natoms = pos.shape[1]
    atom_ppore = natoms // npores

    deviation = np.zeros([nT, npores, atom_ppore])
    for f in tqdm.tqdm(range(nT)):
        for i in range(atom_ppore):
            for j in range(npores):
                deviation[f, j, i] = np.linalg.norm(pos[f, j*atom_ppore + i, :2] - pcenters[f, j, :])

    #deviation = np.reshape(deviation, (nT, natoms))

    fr = np.zeros([nT])
    frstd = np.zeros([nT])
    #
    # for i in range(nT):
    #     fr[i] = np.mean(deviation[i, :])  # + np.std(deviation[i, :]) # maybe?
    #     frstd[i] = np.std(deviation[i, :])

    radii = np.zeros([nT, npores])

    for t in range(nT):
        for p in range(npores):
            radii[t, p] = np.mean(deviation[t, p, :])

    return radii


def put_in_box(pt, x_box, y_box, m, angle):
    """
    :param pt: The point to place back in the box
    :param x_box: length of box in x dimension
    :param y_box: length of box in y dimension
    :param m: slope of box vector
    :param angle: angle between x axis and y box vector
    :return: coordinate shifted into box
    """

    b = - m * x_box  # y intercept of box vector that does not pass through origin (right side of box)

    if pt[1] < 0:
        pt[:2] += [np.cos(angle)*x_box, np.sin(angle)*x_box]  # if the point is under the box
    if pt[1] > y_box:
        pt[:2] -= [np.cos(angle)*x_box, np.sin(angle)*x_box]
    if pt[1] > m*pt[0]:  # if the point is on the left side of the box
        pt[0] += x_box
    if pt[1] < (m*pt[0] + b):  # if the point is on the right side of the box
        pt[0] -= x_box

    return pt


def trace_pores(pos, box, npoints, npores=4, progress=True):
    """
    Find the line which traces through the center of the pores
    :param pos: positions of atoms used to define pore location (args.ref) [natoms, 3]
    :param box: xy box vectors, [2, 2], mdtraj format
    :param npoints: number of points for spline in each pore
    :param npores: number of pores in unit cell (assumed that atoms are number sequentially by pore. i.e. pore 1 atom
    numbers all precede those in pore 2)
    :param progress: set to True if you want a progress bar to be shown

    :return: points which trace the pore center
    """

    single_frame = False
    if np.shape(pos.shape)[0] == 2:
        pos = pos[np.newaxis, ...]  # add a new axis if we are looking at a single frame
        box = box[np.newaxis, ...]
        single_frame = True

    nframes = pos.shape[0]
    atoms_p_pore = int(pos.shape[1] / npores)  # atoms in each pore

    v = np.zeros([nframes, 4, 2])  # vertices of unitcell box
    bounds = []

    v[:, 0, :] = [0, 0]
    v[:, 1, 0] = box[:, 0, 0]
    v[:, 3, :] = np.vstack((box[:, 1, 0], box[:, 1, 1])).T
    v[:, 2, :] = v[:, 3, :] + np.vstack((box[:, 0, 0], np.zeros([nframes]))).T
    center = np.vstack((np.mean(v[..., 0], axis=1), np.mean(v[..., 1], axis=1), np.zeros(nframes))).T

    for t in range(nframes):
        bounds.append(mplPath.Path(v[t, ...]))  # create a path tracing the vertices, v

    angle = np.arcsin(box[:, 1, 1]/box[:, 0, 0])  # specific to case where magnitude of x and y box lengths are equal
    angle = np.where(box[:, 1, 0] < 0, angle + np.pi / 2, angle)  # haven't tested this well yet

    m = (v[:, 3, 1] - v[:, 0, 1]) / (v[:, 3, 0] - v[:, 0, 0])  # slope from points connecting first and third vertices

    centers = np.zeros([nframes, npores, npoints, 3])
    bin_centers = np.zeros([nframes, npores, npoints])

    for t in tqdm.tqdm(range(nframes), disable=(not progress)):
        for p in range(npores):

            pore = pos[t, p*atoms_p_pore:(p+1)*atoms_p_pore, :]  # coordinates for atoms belonging to a single pore

            while np.min(pore[:, 2]) < 0 or np.max(pore[:, 2]) > box[t, 2, 2]:  # because cross-linked configurations can extend very far up and down

                pore[:, 2] = np.where(pore[:, 2] < 0, pore[:, 2] + box[t, 2, 2], pore[:, 2])
                pore[:, 2] = np.where(pore[:, 2] > box[t, 2, 2], pore[:, 2] - box[t, 2, 2], pore[:, 2])

            _, bins = np.histogram(pore[:, 2], bins=npoints)  # bin z-positions

            section_indices = np.digitize(pore[:, 2], bins)  # list that tells which bin each atom belongs to
            bin_centers[t, p, :] = [(bins[i] + bins[i + 1])/2 for i in range(npoints)]

            for l in range(1, npoints + 1):

                atom_indices = np.where(section_indices == l)[0]

                before = pore[atom_indices[0], :]  # choose the first atom as a reference

                shift = transform.translate(pore[atom_indices, :], before, center[t, :])  # shift everything to towards the center

                for i in range(shift.shape[0]):  # check if the points are within the bounds of the unitcell
                    if not bounds[t].contains_point(shift[i, :2]):
                        shift[i, :] = put_in_box(shift[i, :], box[t, 0, 0], box[t, 1, 1], m[t], angle[t])  # if its not in the unitcell, shift it so it is

                c = [np.mean(shift, axis=0)]

                centers[t, p, l - 1, :] = transform.translate(c, center[t, :], before)  # move everything back to where it was

                if not bounds[t].contains_point(centers[t, p, l - 1, :]):  # make sure everything is in the box again
                    centers[t, p, l - 1, :] = put_in_box(centers[t, p, l - 1, :], box[t, 0, 0], box[t, 1, 1], m[t], angle[t])

    if single_frame:
        return centers[0, ...]  # doesn't return bin center yet
    else:
        return centers, bin_centers


def center_of_mass(pos, mass_atoms):
    """ Calculate center of mass of residues over a trajectory

    :param pos: xyz coordinates of atoms
    :param mass_atoms : mass of atoms in order they appear in pos

    :type pos: np.array (nframes, natoms, 3)
    :type mass_atoms: list

    :return: center of mass of each residue at each frame
    """

    nframes = pos.shape[0]
    natoms = len(mass_atoms)

    com = np.zeros([nframes, pos.shape[1] // natoms, 3])  # track the center of mass of each residue

    for f in range(nframes):
        for i in range(com.shape[1]):
            w = (pos[f, i * natoms:(i + 1) * natoms, :].T * mass_atoms).T  # weight each atom in the residue by its mass
            com[f, i, :] = np.sum(w, axis=0) / sum(mass_atoms)  # sum the coordinates and divide by the mass of the residue

    return com


def compdensity(coord, pore_centers, box, cut=1.5, nbins=50, spline=False):
    """ Measure the density of a component as a function of the distance from the pore centers

    :param coord: the coordinates of the component(s) which you want a radial distribution of at each frame
    :param pore_centers: a numpy array of the locations of each pore center at each trajectory frame
    :param cut: cutoff distance for distance calculations. Will not count anything further than cut from the pore center
    :param pores: number of pores (int) default=4
    :param rmax: maximum distance from pore center to calculate density for, default = 3.5 nm
    :param buffer: percentage used to define the location of z planes between which component density will be computed,
           float, default = 0 (i.e. no buffer). Should be between 0 and 1. e.g. for 1 percent, use 0.01 as the buffer

    :type component: numpy.ndarray
    :type pore_centers: numpy.ndarray
    :type cut: float
    :type pores: int
    :type rmax: float
    :type buffer: float

    :return: the density of "component" as a function the distance from the pore center. Also
             returns the calculated bin width for plotting
    """

    nT = coord.shape[0]
    zbox = np.mean(box[:, 2, 2])
    pores = pore_centers.shape[1]
    density = np.zeros([nT, nbins])  # number / nm^3

    if spline:

        npts = pore_centers.shape[2]  # number of points making up the spline in each pore

        edges = np.zeros([nT, pores, npts + 1])  # bin edges, where bin centers are defined by point in spline.
        for t in tqdm.tqdm(range(nT), unit=' Frames'):
            for p in range(pores):
                edges[t, p, 1:-1] = ((pore_centers[t, p, 1:, 2] - pore_centers[t, p, :-1, 2]) / 2) + pore_centers[t, p, :-1, 2]
                edges[t, p, -1] = box[t, 2, 2]

                while np.min(coord[t, :, 2]) < 0 or np.max(coord[t, :, 2]) > box[t, 2, 2]:  # because cross-linked configurations can extend very far up and down
                    coord[t, :, 2] = np.where(coord[t, :, 2] < 0, coord[t, :, 2] + box[t, 2, 2], coord[t, :, 2])
                    coord[t, :, 2] = np.where(coord[t, :, 2] > box[t, 2, 2], coord[t, :, 2] - box[t, 2, 2], coord[t, :, 2])

                zbins = np.digitize(coord[t, :, 2], edges[t, p, :])

                # handle niche case where coordinate lies exactly on the upper or lower bound
                zbins = np.where(zbins == 0, zbins + 1, zbins)
                zbins = np.where(zbins == edges.shape[2], zbins - 1, zbins)

                distances = np.linalg.norm(coord[t, :, :2] - pore_centers[t, p, zbins - 1, :2], axis=1)

                indices = np.where(distances < cut)[0]

                hist, bin_edges = np.histogram(distances[indices], bins=nbins, range=(0, cut))

                density[t, :] += hist

    else:

        for t in tqdm.tqdm(range(nT), unit=' Frames'):
            for p in range(pores):
                # narrow down the positions to those that are within 'cut' of at least one pore
                distances = np.linalg.norm(coord[t, :, :2] - pore_centers[t, p, :], axis=1)
                indices = np.where(distances < cut)[0]
                hist, bin_edges = np.histogram(distances[indices], bins=nbins, range=(0, cut))  # the range option is necessary
                #  to make sure we have equal sized bins on every iteration

                density[t, :] += hist

    # normalize based on volume of anulus where bin is located
    r = np.zeros([nbins])
    for i in range(nbins):
        density[:, i] /= (np.pi * (bin_edges[i + 1] ** 2 - bin_edges[i] ** 2))
        r[i] = (bin_edges[i + 1] + bin_edges[i]) / 2  # center of bins

    density /= (zbox*pores)   # normalize by pore and z-dimension

    return r, density


def minimum_image_distance(dist, box):
    """ Calculate minimum image distances from a vector of distances. This assumes a monoclinic unit cell where the x
    box vector is fixed along the x-axis, the z-box vector is perpendicular to the xy plane, and the y-box vector makes
    an angle, theta, with the x-axis.

    :param d: a vector of distances (n, 3) where n is number of points
    :param box: box vectors meant to enclose d, mdtraj format: (3, 3)

    :return:
    """

    x_box = box[0, 0]  # length of x-box vector
    y_box = box[1, 1]  # perpendicular distance from x-axis to top of box in y-direction
    z_box = box[2, 2]  # length of z-box vector
    d = np.copy(dist)
    angle = np.arcsin(y_box / x_box)  # angle between y-box vector and x-box vector in radians

    # check x coordinates
    while np.max(np.abs(d[:, 0])) > 0.5*x_box:  # iterate in case subtracting/adding box vector length once isn't enough
        d[:, 0] = np.where(d[:, 0] > 0.5*x_box, d[:, 0] - x_box, d[:, 0])
        d[:, 0] = np.where(d[:, 0] < -0.5*x_box, d[:, 0] + x_box, d[:, 0])

    # check y coordinates
    while np.amax(np.abs(d[:, 1])) > 0.5*y_box:  # written differently because np.where didn't know how to handle 2 axes
        d[np.where(d[:, 1] > 0.5*y_box)[0], :2] -= [x_box*np.cos(angle), y_box]
        d[np.where(d[:, 1] < -0.5*y_box)[0], :2] += [x_box*np.cos(angle), y_box]

    # check z coordinates
    while np.max(np.abs(d[:, 2])) > 0.5*z_box:
        d[:, 2] = np.where(d[:, 2] > 0.5*z_box, d[:, 2] - z_box, d[:, 2])
        d[:, 2] = np.where(d[:, 2] < -0.5*z_box, d[:, 2] + z_box, d[:, 2])

    return d


def partition(com, pore_centers, r, buffer=0, unitcell=None, npores=4, spline=False):
    """ Partition residue center of masses into tail and pore region

    :param com: positions of centers of mass of particle whose partition we are calculating
    :param pore_centers: positions of pore centers
    :param r: pore radius, outside of which atoms will be considered in the tail region
    :param buffer: z distance (nm) to cut out from top and bottom of membrane (in cases where there is a water gap)
    :param unitcell: unitcell vectors in mdtraj format (t.unitcell_vectors). Only needed if buffer is used
    :param npores: number of pores
    :param spline: calculate partition with respect to pore spline

    :type com: numpy.ndarray (nT, ncom, 3)
    :type pore_centers: numpy.ndarray (nT, npores, 2) or (nT, npores, 3) or (nT, npores, npts, 3) if spline=True where
    npts=number of points in spline
    :type r: float
    :type buffer: float
    :type unitcell: numpy.ndarray (nT, 3, 3)
    :type npores: int
    :type spline: bool
    """

    nT = com.shape[0]

    if spline:
        npts = pore_centers.shape[2]  # number of points in each spline

    part = np.zeros([nT, com.shape[1]], dtype=bool)  # Will be changed to True if solute in pores

    print('Calculating solute partition...')
    for i in tqdm.tqdm(range(nT)):

        if buffer > 0:
            xy_positions = com[i, (com[i, :, 2] > buffer) & (com[i, :, 2] < unitcell[i, 2, 2] - buffer), :2]
        else:
            xy_positions = com[i, :, :2]

        if spline:

            z = com[i, :, 2]  # extract z-coordinates for this frame
            zbox = unitcell[i, 2, 2]  # z-box vector for this frame

            # make sure z-component of every particle in the box
            while np.max(z) > zbox or np.min(z) < 0:  # might need to do this multiple times
                z = np.where(z > zbox, z - zbox, z)
                z = np.where(z < 0, z + zbox, z)

            zbins = np.digitize(z, np.linspace(0, zbox, npts + 1))

            # handle niche case where coordinate lies exactly on the upper or lower bound
            zbins = np.where(zbins == 0, zbins + 1, zbins)
            zbins = np.where(zbins == npts + 1, zbins - 1, zbins)
            zbins -= 1  # digitize numbers bins starting at 1 (0 is below the bottom bin)

            pore = []

            for p in range(npores):
                d = np.linalg.norm(xy_positions - pore_centers[i, p, zbins, :2], axis=1)
                pore += np.where(d <= r)[0].tolist()

        else:

            pore = []

            for p in range(npores):
                d = np.linalg.norm(xy_positions - pore_centers[i, p, :], axis=1)
                pore += np.where(d <= r)[0].tolist()

        part[i, pore] = True

    return part


def wrap_box(positions, box):
    """ Put all atoms in box

    :param positions: xyz atomic position [n_atoms, 3]
    :param box: box vectors [3, 3] (as obtained from mdtraj t.unitcell_vectors)

    :type positions: np.ndarray
    :type box: np.ndarray

    :return: positions moved into box
    """

    xy = positions[:, :2]  # xy coordinates have dependent changes so this makes things neater below
    z = positions[:, 2]

    xbox, ybox, zbox = box[0, 0], box[1, 1], box[2, 2]

    angle = np.arcsin(ybox / xbox)  # angle between y-box vector and x-box vector in radians
    m = np.tan(angle)
    b = - m * xbox  # y intercept of box vector that does not pass through origin (right side of box)

    while max(xy[:, 1]) > ybox or min(xy[:, 1]) < 0:
        xy[np.where(xy[:, 1] > ybox)[0], :2] -= [xbox*np.cos(angle), ybox]
        xy[np.where(xy[:, 1] < 0)[0], :2] += [xbox * np.cos(angle), ybox]

    while len(np.where(xy[:, 0] < (xy[:, 1] / m))[0]) > 0 or len(np.where(xy[:, 0] > ((xy[:, 1] - b) / m))[0]) > 0:
        xy[np.where(xy[:, 0] < (xy[:, 1] / m))[0], 0] += xbox
        xy[np.where(xy[:, 0] > ((xy[:, 1] - b) / m))[0], 0] -= xbox

    # check z coordinates
    while np.max(z) > zbox or np.min(z) < 0:  # might need to do this multiple times
        z = np.where(z > zbox, z - zbox, z)
        z = np.where(z < 0, z + zbox, z)

    return np.concatenate((xy, z[:, np.newaxis]), axis=1)
