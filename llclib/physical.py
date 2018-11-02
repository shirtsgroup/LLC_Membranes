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


def avg_pore_loc(npores, pos, buffer=0):
    """ Calculate average pore location for each pore at each frame

    :param no_pores: the number of pores in the unit cell
    :param pos: the coordinates of the component(s) which you are using to locate the pore centers
    :param buffer: fraction (of membrane thickness) of top and bottom of membrane to exclude from p2p calculations

    :type no_pores: int
    :type pos: numpy.ndarray, shape(ncomponents, 3) or numpy.ndarray, shape(nframes, ncomponents, 3)
    :type buffer: float

    :return: numpy array containing the x, y coordinates of the center of each pore at each frame
    """

    # Find the average location of the pores w.r.t. x and y

    if len(pos.shape) == 3:  # multiple frames

        nT = np.shape(pos)[0]
        comp_ppore = np.shape(pos)[1] // npores

        #p_center = np.zeros([2, npores, nT])
        p_center = np.zeros([nT, npores, 2])

        for i in range(nT):
            zmax = np.amax(pos[i, :, 2])  # maximum z value for this frame
            zmin = np.amin(pos[i, :, 2])  # minimum z value for this frame
            thick = zmax - zmin
            zmax -= buffer*thick
            zmin += buffer*thick
            for j in range(npores):
                count = 0
                for k in range(comp_ppore*j, comp_ppore*(j + 1)):
                    if zmax >= pos[i, k, 2] >= zmin:
                        p_center[i, j, :] += pos[i, k, :2]
                        count += 1
                p_center[i, j, :] /= count  # take the average

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

    nT = pcenters.shape[2]
    npores = pcenters.shape[1]
    natoms = pos.shape[1]
    atom_ppore = natoms // npores

    deviation = np.zeros([nT, npores, atom_ppore])
    for f in tqdm.tqdm(range(nT)):
        for i in range(atom_ppore):
            for j in range(npores):
                deviation[f, j, i] = np.linalg.norm(pos[f, j*atom_ppore + i, :2] - pcenters[:, j, f])

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
    if pt[1] < m*(pt[0] - b):  # if the point is on the right side of the box
        pt[0] -= x_box

    return pt


def trace_pores(pos, box, npoints, npores=4):
    """
    Find the line which traces through the center of the pores
    :param pos: positions of atoms used to define pore location (args.ref) [natoms, 3]
    :param box: xy box vectors, [2, 2], mdtraj format
    :param npoints: number of points for spline in each pore
    :param npores: number of pores in unit cell (assumed that atoms are number sequentially by pore. i.e. pore 1 atom
    numbers all precede those in pore 2)
    :return: points which trace the pore center
    """

    # assumes trajectory. Probably can add a new axis if a single frame is used
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

    angle = np.arccos(box[:, 1, 1]/box[:, 0, 0])
    angle = np.where(box[:, 1, 0] < 0, angle + np.pi / 2, angle)  # haven't tested this well yet

    m = (v[:, 3, 1] - v[:, 0, 1]) / (v[:, 3, 0] - v[:, 0, 0])  # slope from points connecting first and third vertices

    centers = np.zeros([nframes, npores, npoints, 3])
    bin_centers = np.zeros([nframes, npores, npoints])

    for t in tqdm.tqdm(range(nframes)):
        for p in range(npores):

            pore = pos[t, p*atoms_p_pore:(p+1)*atoms_p_pore, :]  # coordinates for atoms belonging to a single pore
            _, bins = np.histogram(pore[:, 2], bins=npoints)  # bin z-positions
            section_indices = np.digitize(pore[:, 2], bins)  # list that tells which bin each atom belongs to
            bin_centers[t, p, :] = [(bins[i] + bins[i + 1])/2 for i in range(npoints)]

            for l in range(1, npoints + 1):

                atom_indices = np.where(section_indices == l)[0]

                before = pore[atom_indices[0], :]  # choose the first atom as a reference

                shift = transform.translate(pore[atom_indices, :], before, center[t, :])  # shift everything to towards the center

                for i in range(shift.shape[0]):  # check if the points are within the bounds of the unitcell
                    if not bounds[t].contains_point(shift[i, :2]):
                        shift[i, :] = put_in_box(shift[i, :], box[0, 0], box[1, 1], m, angle)  # if its not in the unitcell, shift it so it is

                c = [np.mean(shift, axis=0)]

                centers[t, p, l - 1, :] = transform.translate(c, center[t, :], before)  # move everything back to where it was

                if not bounds[t].contains_point(centers[t, p, l - 1, :]):  # make sure everything is in the box again
                    centers[t, p, l - 1, :] = put_in_box(centers[t, p, l - 1, :], box[0, 0], box[1, 1], m, angle)

    return centers, bin_centers


def center_of_mass(pos, matoms):
    """ Calculate center of mass of residues over a trajectory

    :param pos: xyz coordinates of atoms
    :param names: names of atoms in order they appear in pos
    :param residue: residue object (llclib.topology.Residue())

    :type pos: np.array (nframes, natoms, 3)
    :type names: list
    :type residue: llclib.topology.Residue() object

    :return: center of mass of each residue at each frame
    """

    nframes = pos.shape[0]
    natoms = len(matoms)

    com = np.zeros([nframes, pos.shape[1] // natoms, 3])  # track the center of mass of each residue

    for f in range(nframes):
        for i in range(com.shape[1]):
            w = (pos[f, i * natoms:(i + 1) * natoms, :].T * matoms).T  # weight each atom in the residue by its mass
            com[f, i, :] = np.sum(w, axis=0) / sum(matoms)  # sum the coordinates and divide by the mass of the residue

    return com
