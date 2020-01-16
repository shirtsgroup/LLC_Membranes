#!/usr/bin/env python

import numpy as np
from LLC_Membranes.llclib import file_rw, transform, topology

charge = dict() #same as charge
charge['NA'] = 1.0000

def z_correlation(z, L, v=0.1):
    """ Calculate where to place monomers on the z-axis so that a given correlation length is obtained

    :param z: mean z-positions where monomers will be placed with gaussian probability
    :param L: desired correlation length
    :param v: variance in z position of monomer head groups

    :type z: np.array
    :type L: float
    :type v: float

    :return: locations [np.array[nlayers])
    """

    n = z.shape[0]
    cov = np.zeros([n, n])  # initialize covariance matrix

    decay = v*np.exp(-z / L)  # decay of covariance

    for i in range(z.shape[0]):
        cov[i, i:] = decay[:(n - i)]
        cov[i:, i] = decay[:(n - i)]

    locations = np.random.multivariate_normal(z, cov)

    return locations


class BuildHexagonal:

    def __init__(self, name, npores, p2p, pore_alpha, pore_radius, tilt=0):
        """Initialize geometry of columnar pore structure

        :param name: name of monomer with which the system will be built.
        :param npores: number of pores in the system
        :param p2p: absolute pore-to-pore distance
        :param pore_alpha: angle between x and y box vector. For example if pore_alpha = 120 or 60, you'll get \
        hexagonally packed pores
        :param pore_radius: distance from pore center to place monomer head groups (nm)
        :param tilt: tilt monomer head group with respect to xy plane

        :type name: str
        :type npores: int
        :type p2p: float
        :type pore_alpha: float
        :type pore_radius: float
        :type tilt: float
        """

        # can't use inheritance if there are two LCs
        if type(name) is list:
            self.LC = [topology.LC(i) for i in name]
        else:
            self.LC = [topology.LC(name)]

        self.tilt = tilt
        self.xyz = np.zeros([0, 3])
        self.names = []
        self.all_residues = []
        self.pore_radius = pore_radius

        # currently only implemented for 4 pores
        self.pore_centers = np.zeros([npores, 2])

        pore_alpha_radians = pore_alpha * (np.pi / 180)

        self.pore_centers[1, :] = [p2p*np.cos(pore_alpha_radians), p2p*np.sin(pore_alpha_radians)]
        self.pore_centers[2, :] = [p2p*np.cos(pore_alpha_radians) + p2p, p2p*np.sin(pore_alpha_radians)]
        self.pore_centers[3, :] = [p2p, 0]

        # center pores in box
        self.pore_centers[:, 0] += (p2p / 2) * (1 + np.cos(pore_alpha_radians))
        self.pore_centers[:, 1] += (p2p / 2) * np.sin(pore_alpha_radians)

    def build_column(self, pore, z, theta, correlation=True, var=0, correlation_length=0, pd=0, random_shift=True,
                     mole_fraction=(1.,)):
        """ Place a column at angle theta on xy plane with respect to a pore center

        :param pore: pore number (0 : npores - 1)
        :param z: mean z-positions of monomers in column
        :param theta: angle, with respect to pore center where column should be placed (degrees)
        :param correlation: adjust z positions so there is a correlation length
        :param var: variance in multivariate normal distribution used to make correlated points
        :param correlation_length: length for which correlation between stacked monomers to persist
        :param pd: Angle of wedge created between vertically adjacent monomers. Defined by angle between vectors \
        extending from pore center to monomer head groups.
        :param random_shift: if True, randomly shift columns in z-direction by choosing a displacement from a uniform \
        distribution bounded by (0, d), where d is the vertical distance between stacked monomers
        :param seed: random seed if you want to reproduce randomly displaced structures
        :param mole_fraction: mol fraction of each type of monomer. This only has meaning if the system is built with a
        mixture of monomers.

        :type pore: int
        :type z: np.array
        :type theta: float
        :type correlation: bool
        :type var: float
        :type correlation_length: float
        :type pd: float
        :type random_shift: bool
        :type mole_fraction: tuple of floats or numpy.ndarray
        """

        if correlation:
            z = z_correlation(z, correlation_length, v=var)
            if random_shift:
                dbwl = z[1] - z[0]  # distance between stacked monomers
                z += np.random.uniform(0, dbwl)

        elif random_shift:
            dbwl = z[1] - z[0]  # distance between stacked monomers
            z += np.random.uniform(0, dbwl)

        displaced_theta = pd

        # natoms = self.LC_positions.shape[0]  # number of atoms including ions
        # pos = np.copy(self.LC_positions)
        # pos[:, 0] += self.pore_radius
        # displaced_pos = np.copy(pos)

        natoms = [lc.LC_positions.shape[0] for lc in self.LC]
        pos = []
        for i, lc in enumerate(self.LC):
            pos.append(np.copy(lc.LC_positions))
            pos[i][:, 0] += self.pore_radius
        displaced_pos = np.copy(pos)

        pos = [transform.rotate_coords_z(p, theta) for p in pos]
        displaced_pos = [transform.rotate_coords_z(p, theta + displaced_theta) for p in displaced_pos]

        if len(mole_fraction) != len(self.LC):
            raise Exception('You must supply one mole fraction for each monomer.')

        mole_fraction = np.array(mole_fraction)
        mole_fraction /= mole_fraction.sum()

        monomer_ids = np.random.choice(np.arange(len(self.LC)), p=mole_fraction, size=z.size)

        atom_indices = np.cumsum([0] + [natoms[i] for i in monomer_ids])
        column = np.zeros([atom_indices[-1], 3])

        before = np.array([0, 0, 0])
        for l in range(z.size):
            mon = monomer_ids[l]
            if l % 2 == 0:
                column[atom_indices[l]: atom_indices[l + 1]] = transform.translate(pos[mon], before, np.array([0, 0, z[l]]))
            else:
                column[atom_indices[l]: atom_indices[l + 1]] = transform.translate(displaced_pos[mon], before,
                                                                             np.array([0, 0, z[l]]))
            self.names += self.LC[mon].LC_names
            self.all_residues += self.LC[mon].LC_residues

        column = transform.translate(column, before, np.array([self.pore_centers[pore, 0], self.pore_centers[pore, 1], 0]))

        self.xyz = np.concatenate((self.xyz, column))

    def write_gro(self, out, ucell):
        """ Write coordinate file in .gro format

        :param out: name of output .gro file
        :param ucell: unitcell vectors

        :type out: str
        :type ucell: np.ndarray, shape(3,3)

        """

        file_rw.write_gro_pos(self.xyz, out, ids=self.names, res=self.all_residues, ucell=ucell)

    def reorient_monomer(self):
        """ Align monomer head groups so they are coplanar with the xy plane and vector pointing towards pore center is
        aligned with x-axis
        """

        for lc in self.LC:
            lc.align_monomer(tilt=self.tilt)

    # def align_plane(self):
    #     """ Align the atoms defined by the plane_indices attribute of LC with the xy plane """
    #
    #     plane_atoms = np.zeros([3, 3])
    #     for i in range(plane_atoms.shape[0]):
    #         plane_atoms[i, :] = self.LC_positions[self.plane_indices[i], :]
    #
    #     R = transform.rotateplane(plane_atoms, angle=self.tilt)  # generate rotation matrix
    #
    #     b = np.ones([1])
    #     for i in range(self.LC_positions.shape[0]):
    #         coord = np.concatenate((self.LC_positions[i, :], b))
    #         x = np.dot(R, coord)
    #         self.LC_positions[i, :] = x[:3]
    #
    # def translate_to_origin(self):
    #     """ Translate molecule to the origin using the ref_atom_index attribute of LC """
    #
    #     self.LC_positions = transform.translate(self.LC_positions, self.LC_positions[self.ref_atom_index, :], [0, 0, 0])
    #
    # def align_with_x(self):
    #     """ Align vector defined by lineatoms in LC object with x axis """
    #
    #     v = np.array([self.LC_positions[self.lineatoms[0], :2] - self.LC_positions[self.lineatoms[1], :2]])
    #     angle = np.arctan2(v[0, 0, 1], v[0, 0, 0])
    #     self.LC_positions = transform.rotate_coords_z(self.LC_positions, - angle * 180 / np.pi)

    def reorder(self):
        """ reorder coordinate, residues and atom names so that all residues are separated. Some scripts (might) still
        rely on this ordering.
        """

        residues = []
        for lc in self.LC:
            for r in lc.residues:
                residues.append(r)

        ordered = []
        for r in residues:
            ndx = [i for i, a in enumerate(self.all_residues) if a == r]
            ordered += ndx
        #     residue_indices.append(ndx)
        #
        # ordered = []
        # for i in range(len(residue_indices)):
        #     ordered += residue_indices[i]

        self.xyz = self.xyz[ordered, :]
        self.all_residues = [self.all_residues[i] for i in ordered]
        self.names = [self.names[i] for i in ordered]

        def exchange_ion(self, original_ion, new_ion):
        """ This function will switch out the ion from the system and replace with a new ion.
"""
        #modify self.xyz (coordinates), self.names, self.all_residues get rid of all NA stuff and replace with something else.

        charge[new_ion] #figuring out how many ions needed
	#use some sodium position charges... for +2, use every second one. for +3, every third NA coordinate

	#how to manipulate numpy array to do what I think I want it to do.
