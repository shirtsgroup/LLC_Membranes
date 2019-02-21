#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
from LLC_Membranes.llclib import topology, physical, file_rw
from LLC_Membranes.analysis import solute_partitioning
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Measure the nematic order parameter based on the vectors extending'
                                                 'from a solute COM to the pore center and in the direction of the'
                                                 'of the solute as defined in the coordinate file.')

    parser.add_argument('-t', '--traj', default='PR.xtc', help='Name of system trajectory file (.trr or .xtc)')
    parser.add_argument('-g', '--gro', default='PR.gro', help='Name of structure file where water will be added')
    parser.add_argument('-b', '--build_monomer', default='NAcarb11V.gro', help='Name of monomer structure file used to'
                                                                               ' build LLC membrane')
    parser.add_argument('-s', '--solute', default='ETH', help='Name of solute to be analysed')
    parser.add_argument('-l', '--load', default=False, help='Load pickled object')
    parser.add_argument('-w', '--wall_location', default=False, type=float, help='Distance from pore center where '
                                                                                 'solute center of mass is considered '
                                                                                 'in the head group region')

    return parser

# this works but the variable names are less clear
# class System(topology.Solute):
#
#     def __init__(self, solute, gro, traj, build_monomer, spline=False):
#
#         super().__init__(solute)


def moving_average(a, n=3):

    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]

    return ret[n - 1:] / n


class System(object):

    def __init__(self, solute, gro, traj, build_monomer, spline=False, npores=4):

        self.npores = npores
        self.traj = traj
        self.gro = gro
        self.build_monomer = build_monomer

        print('Loading trajectory...', end='', flush=True)
        self.t = md.load(self.traj, top=self.gro)
        print('Done!')

        self.solute = topology.Solute(solute)  # Create solute object defining solute direction vector
        self.monomer = topology.LC(self.build_monomer)  # Object for liquid crystal monomer used to build system

        # monomer head group only for testing purposes
        # ref = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
        # solutes_indices = [a.index for a in self.t.topology.atoms if a.residue.name == solute and a.name in ref]
        # self.nsolute = 400
        # self.solute_vectors = self.direction_vectors()
        # self.com = physical.center_of_mass(self.t.xyz[:, solutes_indices, :], [v for v in self.solute.mass.values()][:6])

        solutes_indices = [a.index for a in self.t.topology.atoms if a.residue.name == solute]  # all indices of solute
        self.nsolute = len(solutes_indices) // self.solute.natoms  # number of solutes in system

        self.com = None
        self.solute_vectors = None
        self.pore_centers = None
        self.pore_vectors = None
        self.costheta = None
        self.nematic_order_parameter = None
        self.solute_partition = None

    def direction_vectors(self):

        print('Calculating direction vector of %s using %s and %s' % (self.solute.name,
                                                                      self.solute.direction_vector[0][0],
                                                                      self.solute.direction_vector[1][0]))

        self.solute_vectors = np.zeros([self.t.n_frames, self.nsolute, 2])
        self.com = np.zeros([self.t.n_frames, self.nsolute, 3])
        # self.com = physical.center_of_mass(self.t.xyz[:, solutes_indices, :], [v for v in self.solute.mass.values()])

        back_atom = [a.index for a in self.t.topology.atoms if a.residue.name == self.solute.name and
                      a.name in self.solute.direction_vector[0]]

        front_atom = [a.index for a in self.t.topology.atoms if a.residue.name == self.solute.name and
                      a.name in self.solute.direction_vector[1]]

        for t in range(self.t.n_frames):
            self.solute_vectors[t, ...] = self.t.xyz[t, front_atom, :2] - self.t.xyz[t, back_atom, :2]  # xy projection only
            self.com[t, ...] = (self.t.xyz[t, front_atom, :] + self.t.xyz[t, back_atom, :]) / 2

    def calculate_pore_centers(self, spline=True):

        reference_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms and
                           a.residue.name in self.monomer.residues]  # atoms whose avg positions to define pore centers

        self.pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, reference_atoms, :], self.t.unitcell_vectors,
                                                  spline=spline, progress=True)  # get spline running through pores

    def pore_center_vectors(self):

        self.pore_vectors = np.zeros([self.t.n_frames, self.nsolute, 2])

        solute_per_pore = self.nsolute // self.npores

        for t in range(self.t.n_frames):
            for p in range(self.npores):
                v_centers = physical.wrap_box(self.com[t, p*solute_per_pore:(p+1)*solute_per_pore, :],
                                              self.t.unitcell_vectors[t, ...])  # solute vector centers wrapped into box
                bins = np.digitize(v_centers[:, 2], self.pore_centers[t, p, :, 2]) - 1
                pore_centers = self.interpolate_spline(v_centers[:, 2], bins, self.pore_centers[t, p, ...],
                                                       self.t.unitcell_vectors[t, 2, 2])

                # self.pore_vectors[t, (p*solute_per_pore):(p+1)*solute_per_pore, :] = self.pore_centers[t, p, :] - \
                #                                              self.com[t, (p*solute_per_pore):(p+1)*solute_per_pore, :2]

                # self.pore_vectors[t, (p*solute_per_pore):(p+1)*solute_per_pore, :] = self.pore_centers[t, p, bins, :2] \
                #                                                                      - v_centers[:, :2]
                self.pore_vectors[t, (p*solute_per_pore):(p+1)*solute_per_pore, :] = pore_centers - v_centers[:, :2]

    def interpolate_spline(self, positions, bins, spline, zbox):
        """ Linearly interpolate between points in pore spline to determine the pore center as a function of z

        :param positions:
        :param bins:
        :param spline:
        :param zbox:
        :return:
        """

        bounds = np.concatenate((bins[:, np.newaxis], (bins + 1)[:, np.newaxis]), axis=1)

        # modify spline so that periodic spline copies are included for upper and lower bounds
        spline = np.concatenate((spline, spline[np.newaxis, 0, :], spline[np.newaxis, -1, :]), axis=0)
        spline[-2, 2] += zbox
        spline[-1, 2] -= zbox

        xy = np.zeros_like(bounds, dtype=float)

        lower, upper = spline[bounds[:, 0], :], spline[bounds[:, 1], :]
        bin_width = upper - lower
        zd = positions - lower[:, 2]
        percentage = zd / bin_width[:, 2]

        for i in range(percentage.size):
            xy[i, :] = lower[i, :2] + percentage[i] * bin_width[i, :2]

        # for i, b in enumerate(bounds):
        #     lower, upper = spline[b[0], :], spline[b[1], :]
        #     bin_width = upper - lower
        #     zd = positions[i] - lower[2]
        #     percentage = zd / bin_width[2]
        #     xy[i, :] = lower[:2] + percentage*bin_width[:2]

        return xy

    def angles(self):
        """
        Calculate angle between a vector and plane
        :param vector: vector(s) whose angle with 'plane' we want
        :param plane: vector normal to plane
        :return:
        """

        nT = self.solute_vectors.shape[0]
        self.costheta = np.zeros([nT, self.nsolute])

        for i in range(nT):
            mag_vsolute = np.linalg.norm(self.solute_vectors[i, ...], axis=1)  # magnitude of each solute vector
            mag_vpore = np.linalg.norm(self.pore_vectors[i, ...], axis=1)
            for s in range(self.nsolute):  # since there is no vectorwise dot product apparently
                self.costheta[i, s] = np.dot(self.solute_vectors[i, s, :], self.pore_vectors[i, s, :]) / \
                                 (mag_vsolute[s] * mag_vpore[s])

    def nematic(self):
        """
        Calculate the nematic order parameter of a liquid crystal system
        See p. 168 of : http://www.dsf.unica.it/~fiore/libricorsoptr/Chaikin%20Lubensky%20-%20Principles%20of%20Condensed%20Matter%20Physics.jb2.pdf
        :param costheta: cosine of the angle made with director vector (get from "angles" function above)
        :return: nematic order parameter at each frame
        """

        self.nematic_order_parameter = np.zeros([self.costheta.shape[0]])  # calculate S for each frame
        for i in range(self.costheta.shape[0]):
            #S[i] = np.sum((3 * self.costheta[i, :-2] ** 2 - 1) / 2) # 3D
            if self.solute_partition is not None:
                ndx = np.where(~self.solute_partition[i, ])[0]
            else:
                ndx = np.arange(self.costheta.shape[1])
            # 2D (https://physics.stackexchange.com/questions/65358/2-d-orientational-order-parameter)
            self.nematic_order_parameter[i] = np.sum(2 * self.costheta[i, ndx] ** 2 - 1) / len(ndx)

        return self.nematic_order_parameter

    def plot(self, show=True, save=False, ma=False, ybounds=[-1, 1]):

        fig, ax = plt.subplots(1, 1)

        if ma:
            ax.plot(self.t.time[:-(ma - 1)] / 1000, moving_average(self.nematic_order_parameter, n=ma), linewidth=2)
        else:
            ax.plot(self.t.time / 1000, self.nematic_order_parameter, linewidth=2)

        ax.set_ylabel('Order Parameter', fontsize=14)
        ax.set_xlabel('Time (ns)', fontsize=14)
        ax.set_ylim(ybounds)
        ax.tick_params(axis='both', labelsize=14)
        plt.tight_layout()

        if show:
            plt.show()

    def write_com_positions(self, frame=0, name='com.gro'):

        file_rw.write_gro_pos(self.com[frame, ...], name)

    def write_spline_positions(self, frame=0, name='spline.gro'):

        file_rw.write_gro_pos(self.pore_centers[frame, ...].reshape(self.pore_centers.shape[1] *
                                                                    self.pore_centers.shape[2], 3), name, name='K')

    def partition(self, r):

        self.solute_partition = physical.partition(self.com, self.pore_centers, r, unitcell=self.t.unitcell_vectors,
                                                   spline=True)


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.load:
        sys = file_rw.load_object(args.load)
    else:
        sys = System(args.solute, args.gro, args.traj, args.build_monomer)
        sys.direction_vectors()
        sys.calculate_pore_centers()
        sys.pore_center_vectors()
        if args.wall_location:
            sys.partition(args.wall_location)
        sys.angles()
        sys.nematic()
        file_rw.save_object(sys, 'order.pl')

    sys.nematic()
    sys.plot(ma=False)
