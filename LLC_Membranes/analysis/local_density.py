#!/usr/bin/env python

import mdtraj as md
import numpy as np
from LLC_Membranes.llclib import file_rw
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
import tqdm
import time


class Cube(object):

    def __init__(self, xrange, yrange, zrange):

        self.xrange = xrange  # (xmin, xmax)
        self.yrange = yrange
        self.zrange = zrange

    def contains_point(self, p):

        return all([self.xrange[0] <= p[0] <= self.xrange[1],
                    self.yrange[0] <= p[1] <= self.yrange[1],
                    self.zrange[0] <= p[2] <= self.zrange[1]])


class LocalDensity:

    def __init__(self, traj, gro, start_frame=0):
        """ Find the local density in a monoclinic unit cell
        """

        print('Loading Trajectory...', end='', flush=True)
        t = md.load(traj, top=gro)[start_frame:]
        print('Done!')

        keep = [a.index for a in t.topology.atoms if a.element.name is not 'hydrogen']  # don't include hydrogens

        self.pos = t.xyz[:, keep, :]

        self.box_vectors = t.unitcell_vectors
        box = t.unitcell_vectors.mean(axis=0)
        self.xbox, self.ybox, self.zbox = box[0, 0], box[1, 1], box[2, 2]
        self.xshift = np.sqrt(self.xbox ** 2 - self.ybox ** 2)
        self.cube = None

        self.xrange = None
        self.yrange = None
        self.zrange = None

        self.density = []

        self.n_frames = t.n_frames

    def calculate(self, buffer=1, bins=(10, 10, 10), create_interpolator=True, shrink_box=False):
        """
        :params buffer: minimum distance between edge of monoclinic box and cubic box that'll surround it (nm)
        :params bins: number of bins in each dimension. If int is passed, all dimensions will have the same number of \
        bins.
        :param shrink_box: if True, only atoms within the self.cube object will be kept. This can help with memory, \
        but is generally not necessary

        :type buffer: float
        :type bins: int or tuple
        :type shrink_box: bool
        """

        self.xrange = (-buffer, buffer + self.xshift + self.xbox)
        self.yrange = (-buffer, buffer + self.ybox)
        self.zrange = (-buffer, buffer + self.zbox)

        self.cube = Cube(self.xrange, self.yrange, self.zrange)

        if not create_interpolator:
            self.density = np.zeros([self.n_frames, bins[0], bins[1], bins[2]])

        for t in tqdm.tqdm(range(self.n_frames), disable=False):

            pos = self._replicate_periodically(self.pos[t, ...], shrink_box=shrink_box)

            d, edges = np.histogramdd(pos, bins=bins, range=(self.xrange, self.yrange, self.zrange), density=True)

            if create_interpolator:
                self.density.append(self._interpolator(d, edges))
            else:
                self.density[t, ...] = d

    def _check_points(self, pos):

        keep = []
        for i, p in enumerate(pos):
            if self.cube.contains_point(p):
                keep.append(i)

        return keep

    def _replicate_periodically(self, initial_positions, shrink_box=False):
        """ Replicate positions periodically in +/- x, y and z direction. But only keep points that are in bounding box.

        :param initial_positions: (n, 3) array of positions to be periodically replicated in all directions
        :param shrink_box: if True, only atoms within the self.cube object will be kept. This can help with memory, but
        is generally not necessary

        :type initial_positions: np.ndarray
        :type shrink_box: bool
        """

        shifts = [-1, 0, 1]
        pos = np.empty([0, 3])
        for x in shifts:
            positions = initial_positions
            positionsx = np.copy(positions)
            positionsx[:, 0] += self.xbox * x
            for y in shifts:
                positionsy = np.copy(positionsx)
                positionsy[:, 1] += self.ybox * y
                positionsy[:, 0] += self.xshift * y
                for z in shifts:
                    positionsz = np.copy(positionsy)
                    positionsz[:, 2] += self.zbox * z
                    if shrink_box:
                        keep = self._check_points(positionsz)
                        pos = np.concatenate((pos, positionsz[keep, :]))
                    else:
                        pos = np.concatenate((pos, positionsz))

        return pos

    def _interpolator(self, h, edges, type='regular'):
        """ Use a regular grid interpolator to make an object that allows interpolation at arbitrary points within the
        density histogram.

        :param h: 3D histogram created from np.histogramdd
        :param edges: The edges of the histogram bins. The centers of these bins will be used as the points defining \
        the regular grid
        :param type: type of interpolation to perform ('linear': LinearNDInterpolator, 'regular': RegularGridInterpolator)

        :type h: np.array
        :type edges: list of arrays in each dimension
        :type type: str
        """

        xedge_shift = (edges[0][1] - edges[0][0]) / 2
        yedge_shift = (edges[1][1] - edges[1][0]) / 2
        zedge_shift = (edges[2][1] - edges[2][0]) / 2

        x = [i + xedge_shift for i in edges[0][:-1]]
        y = [i + yedge_shift for i in edges[1][:-1]]
        z = [i + zedge_shift for i in edges[2][:-1]]

        if type == 'regular':
            return RegularGridInterpolator((x, y, z), h)
        elif type == 'linear':
            return LinearNDInterpolator((x, y, z), h)
        else:
            raise Exception('%s interpolation is not implemented' % type)
