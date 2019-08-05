#! /usr/bin/env python

from __future__ import division
from builtins import range
from LLC_Membranes.llclib import file_rw, transform, topology
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
                print(i)

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


def trace_pores(pos, box, npoints, npores=4, progress=True, save=True, savename='spline.pl'):
    """
    Find the line which traces through the center of the pores
    :param pos: positions of atoms used to define pore location (args.ref) [natoms, 3]
    :param box: xy box vectors, [2, 2], mdtraj format
    :param npoints: number of points for spline in each pore
    :param npores: number of pores in unit cell (assumed that atoms are number sequentially by pore. i.e. pore 1 atom
    numbers all precede those in pore 2)
    :param progress: set to True if you want a progress bar to be shown
    :param save: save spline as pickled object

    :return: points which trace the pore center
    """

    try:
        print('Attempting to load spline ... ', end='', flush=True)
        spline = file_rw.load_object(savename)
        print('Success!')

        return spline[0], spline[1]

    except FileNotFoundError:

        print('%s not found ... Calculating spline' % savename)

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
                        while not bounds[t].contains_point(shift[i, :2]):
                            shift[i, :] = put_in_box(shift[i, :], box[t, 0, 0], box[t, 1, 1], m[t], angle[t])  # if its not in the unitcell, shift it so it is

                    c = [np.mean(shift, axis=0)]

                    centers[t, p, l - 1, :] = transform.translate(c, center[t, :], before)  # move everything back to where it was

                    while not bounds[t].contains_point(centers[t, p, l - 1, :]):  # make sure everything is in the box again
                        centers[t, p, l - 1, :] = put_in_box(centers[t, p, l - 1, :], box[t, 0, 0], box[t, 1, 1], m[t], angle[t])

        if single_frame:
            return centers[0, ...]  # doesn't return bin center yet

        else:

            if save:
                file_rw.save_object((centers, bin_centers), savename)

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


def residue_center_of_mass(t, res):
    """ Calculate the center of mass versus time of a residue in an MD trajectory

    :param t: mdtraj trajectory object
    :param res: name of residue to track

    :type t: object
    :type res: str

    :return: center of mass of residue versus time
    """

    residue = topology.Residue(res)  # get resiude attributes

    ndx = [a.index for a in t.topology.atoms if a.residue.name == res]  # index of all residue atoms
    names = [a.name for a in t.topology.atoms if a.residue.name == res][:residue.natoms]  # names of atoms in one residue
    mass = [residue.mass[x] for x in names]  # mass of atoms in order that they appear in file
    print('Calculating center of mass trajectories of residue %s' % residue.name)

    return center_of_mass(t.xyz[:, ndx, :], mass)  # determine center of mass trajectories


def compdensity(coord, pore_centers, box, cut=1.5, nbins=50, spline=False):
    """ Measure the density of a component as a function of the distance from the pore centers.

    :param coord: the coordinates of the component(s) which you want a radial distribution of at each frame
    :param pore_centers: a numpy array of the locations of each pore center at each trajectory frame
    :param cut: cutoff distance for distance calculations. Will not count anything further than cut from the pore center
    :param nbins: number of bins in r direction
    :param spline: calculate RDF with respect to spline

    :type coord: numpy.ndarray
    :type pore_centers: numpy.ndarray
    :type cut: float
    :type nbins: int
    :type spline: bool


    :return: Radial distance from pore center r, and the density of a species, whose positions are defined by
    `coordinates`, as a function the distance from the pore center.
    """

    nT = coord.shape[0]
    pores = pore_centers.shape[1]
    density = np.zeros([nT, nbins])  # number / nm^3

    for t in tqdm.tqdm(range(nT), unit=' Frames'):
        for p in range(pores):

            if spline:
                distances = radial_distance_spline(pore_centers[t, p, ...], coord[t, ...], box[t, ...])
            else:
                distances = np.linalg.norm(coord[t, :, :2] - pore_centers[t, p, :], axis=1)

            hist, bin_edges = np.histogram(distances, bins=nbins, range=(0, cut))

            density[t, :] += hist

        density[t, :] /= (pores * box[t, 2, 2])  # normalize by z-dimension

    # normalize based on volume of anulus where bin is located (just need to divide by area since height done above)
    r = np.zeros([nbins])
    for i in range(nbins):
        density[:, i] /= (np.pi * (bin_edges[i + 1] ** 2 - bin_edges[i] ** 2))
        r[i] = (bin_edges[i + 1] + bin_edges[i]) / 2  # center of bins

    return r, density


def distance_from_pore_center(coord, pore_centers, box, spline=False):
    """ Measure the density of a component as a function of the distance from the pore centers.

    :param coord: the coordinates of the component(s) which you want a radial distribution of at each frame
    :param pore_centers: a numpy array of the locations of each pore center at each trajectory frame
    :param cut: cutoff distance for distance calculations. Will not count anything further than cut from the pore center
    :param

    :type coord: numpy.ndarray
    :type pore_centers: numpy.ndarray
    :type cut: float

    :return: Radial distance of each individual solute/component, defined by coords, as a function of time
    """

    nT = coord.shape[0]
    pores = pore_centers.shape[1]
    nsolute = coord.shape[1]

    r_distances = np.zeros([nT, nsolute])

    for t in tqdm.tqdm(range(nT), unit=' Frames'):
        rd = np.zeros([nsolute, pores])
        for p in range(pores):

            if spline:
                rd[:, p] = radial_distance_spline(pore_centers[t, p, ...], coord[t, ...], box[t, ...])
            else:
                rd[:, p] = np.linalg.norm(coord[t, :, :2] - pore_centers[t, p, :], axis=1)

        # Move the minimum solute--pore-center distance for each solute to the first index of rd
        # This removes any assumption that there is a constant number of solutes per pore and that the solute
        # stays in the same pore.
        for i, r in enumerate(rd):  # there is probably a vectorized way to do this with argsort
            rd[i, :] = r[np.argsort(r)]

        r_distances[t, :] = rd[:, 0]

    return r_distances


def radial_distance_spline(spline, com, box):
    """ Calculate radial distance from pore center based on distance from center of mass to closest z point in spline

    :param spline: coordinates of spline for a single pore and frame
    :param com: atomic center of mass z-coordinates
    :param zbox: z box dimension (nm)

    :type spline: np.ndarray [npts_spline, 3]
    :type com: np.ndarray [n_com, 3]
    :type zbox: float

    :return: array of distances from pore center
    """

    edges = np.zeros([spline.shape[0] + 1])
    edges[1:-1] = ((spline[1:, 2] - spline[:-1, 2]) / 2) + spline[:-1, 2]
    edges[-1] = box[2, 2]

    com = wrap_box(com, box)
    # while np.min(com[:, 2]) < 0 or np.max(com[:, 2]) > zbox:  # because cross-linked configurations can extend very far up and down
    #     com[:, 2] = np.where(com[:, 2] < 0, com[:, 2] + zbox, com[:, 2])
    #     com[:, 2] = np.where(com[:, 2] > zbox, com[:, 2] - zbox, com[:, 2])

    zbins = np.digitize(com[:, 2], edges)

    # handle niche case where coordinate lies exactly on the upper or lower bound
    zbins = np.where(zbins == 0, zbins + 1, zbins)
    zbins = np.where(zbins == edges.size, zbins - 1, zbins)

    return np.linalg.norm(com[:, :2] - spline[zbins - 1, :2], axis=1)


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
    :param unitcell: unitcell vectors in mdtraj format (t.unitcell_vectors). Only needed if buffer and/or spline is used
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

    :return part: boolean numpy array with shape (nT, com.shape[1]) where True indicates a center of mass that is
    inside the inner region (i.e. < r)
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


def wrap_box(positions, box, tol=1e-6):
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

    # added tolerance for corner case
    while len(np.where(xy[:, 0] - (xy[:, 1] / m) < -tol)[0]) > 0 or \
            len(np.where(xy[:, 0] - ((xy[:, 1] - b) / m) > 0)[0]) > 0:
        xy[np.where(xy[:, 0] < (xy[:, 1] / m))[0], 0] += xbox
        xy[np.where(xy[:, 0] > ((xy[:, 1] - b) / m))[0], 0] -= xbox

    # check z coordinates
    while np.max(z) > zbox or np.min(z) < 0:  # might need to do this multiple times
        z = np.where(z > zbox, z - zbox, z)
        z = np.where(z < 0, z + zbox, z)

    return np.concatenate((xy, z[:, np.newaxis]), axis=1)
