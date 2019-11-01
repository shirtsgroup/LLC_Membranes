#!/usr/bin/env python

import numpy as np
from LLC_Membranes.setup import surfaces
from LLC_Membranes.llclib import topology, transform, file_rw
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy import spatial
import random

mw = {'glycerol': 92.09382}


class UnitCellError(Exception):
    """ Raised if invalid phase specified """

    def __init__(self, message):

        super().__init__(message)


class CurvatureError(Exception):
    """ Raised if invalid value of curvature is passed """

    def __init__(self, message):

        super().__init__(message)


class BicontinuousCubicBuild(topology.LC):

    def __init__(self, monomer, space_group, dimensions, weight_percent, density):
        """ Initialize a build of a bicontinuous cubic phase unit cell

        :param monomer: name of monomer with which to build phase
        :param space_group: name of space group into which monomers are arranged (i.e. gyroid, schwarzD etc.)
        :param dimensions: length of each edge of the unit cell (nm). The box is cubic, so they must all be the same \
        (or just one length specified)
        :param weight_percent: percent by weight of monomer in membrane
        :param density: experimentally derived density of membrane

        :type monomer: str
        :type space_group: str
        :type dimensions: float, list of floats
        :type weight_percent: int, float
        :type density: float
        """

        super().__init__(monomer)

        self.space_group = space_group

        if type(dimensions) is list:
            if len(set(dimensions)) > 1:
                raise UnitCellError('Your x, y and z box vectors are not all the same. The unit cell must be a cubic '
                                    'box')
            self.period = dimensions[0]
        else:
            self.period = dimensions
        #
        # print(self.period)
        # exit()

        # calculate nubmer of monomers that will go in the unit cell
        avogadros_number = 6.022 * 10 ** 23
        mass = self.MW / avogadros_number  # mass of a single monomer (g)
        volume = self.period ** 3 * 10 ** -21  # volume of unit cell (cm ^ 3)
        monomer_mass = (weight_percent / 100.) * density * volume  # mass of all monomers in unit cell
        self.nmon = int(monomer_mass / mass)

        # initialize other variables
        self.grid = None
        self.final_positions = None
        self.all_residues = None
        self.all_names = None

    def gen_grid(self, n, curvature, plot=False):
        """ Create an n x n x n grid of points which lie close to the surface describing the unit cell space group

        :param n: number of points in the x, y and z direction
        :param curvature: determines mean curvature of hydrophilic/hydrophobic interface and thus whether the phase is
        normal or inverted. {> 0 : QI phase, < 0: QII phase}'
        :param plot: show a 3D scatter plot of the surface

        :type n: int
        :type curvature: float
        :type plot: bool
        """

        self.grid = surfaces.gridgen(self.space_group, 0, self.period, n, c=curvature)

        if plot:

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(self.grid[:, 0], self.grid[:, 1], self.grid[:, 2])
            plt.show()

    def determine_monomer_placement(self, r=0.4):
        """ Choose grid points where monomers will be placed

        :param r: Any points within this distance of a chosen grid point will be removed as potential new grid points \
        in an effort to keep point sufficiently spaced.

        :type r: float
        """

        count = 0
        new_grid = np.zeros([self.nmon, 3])  # grid of final monomer placement points

        while count < self.nmon:

            delete = random.randint(0, self.grid.shape[0] - 1)  # pick a random point on the grid
            new_grid[count, :] = self.grid[delete, :]  # assign point to new grid

            # find all points within a distance r of the chosen point
            tree = spatial.cKDTree(self.grid)
            nn = tree.query_ball_point(self.grid[delete, :], r)  # all nearest neighbors to grid[delete, :] within radius r
            nn.append(delete)
            self.grid = np.delete(self.grid, nn, 0)  # delete nearest neighbors and itself from grid

            count += 1

            if self.grid.shape[0] == 0:
                break

        if count < self.nmon:
            # probably better to raise an exception here and quit
            print('Only %s monomers were placed. Try again with a lower r value or different number of grid points' % count)
            self.nmon = count

        self.grid = new_grid

    def place_monomers(self, shift=0):
        """ Place monomers perpendicular to the space group surface at grid point locations

        :param shift: translate monomer along vector perpendicular to space group surface by this amount (nm). This \
        parameter effectively controls the pore size

        :type shift: float
        """

        # if curvature not in [-1, 1]:
        #     raise CurvatureError('The value for curvature must be either -1 or 1')

        self.final_positions = np.zeros([self.natoms * self.nmon, 3])

        # multiple atoms can be used to specify the lineatoms and reference atoms. Their average positions are used

        l1 = self.LC_positions[self.lineatoms[1], :].mean(axis=0)
        l2 = self.LC_positions[self.lineatoms[0], :].mean(axis=0)

        linevector = l1 - l2
        reference_position = self.LC_positions[self.ref_atom_index, :].mean(axis=0)

        for i in range(self.nmon):

            n = surfaces.gradient(self.grid[i, :], self.space_group, period=self.period)  # vector normal to surface at point grid[i, :]

            R = transform.Rvect2vect(linevector, n)  # rotation matrix to rotate monomer in same direction as n

            # translate to origin
            xyz_origin = transform.translate(self.LC_positions, reference_position, np.array([0, 0, 0]))

            xyz_origin = transform.rotate_coords(xyz_origin, R)  # rotate all points in bcc monomer with rotation matrix

            ref = xyz_origin[self.ref_atom_index, :].mean(axis=0)  # avg loc of reference atoms changes after rotation

            placement = self.grid[i, :] - shift * n

            xyz_origin = transform.translate(xyz_origin, ref,
                                             placement)  # move monomer to grid point w.r.t. reference point on monomer

            self.final_positions[i * self.natoms:(i + 1) * self.natoms, :] = xyz_origin

    def reorder(self):
        """ reorder coordinate, residues and atom names so that residues are separated """

        ordered = []
        for r in self.residues:
            ordered += [i for i, a in enumerate(self.LC_residues * self.nmon) if a == r]

        self.final_positions = self.final_positions[ordered, :]
        all_residues = self.LC_residues * self.nmon
        all_names = self.LC_names * self.nmon
        self.all_residues = [all_residues[i] for i in ordered]
        self.all_names = [all_names[i] for i in ordered]

    def write_final_configuration(self, name='initial.gro', box=None):

        if not box:
            box = [self.period, self.period, self.period]

        file_rw.write_gro_pos(self.final_positions, name, res=self.all_residues, ids=self.all_names,
                              box=box)
