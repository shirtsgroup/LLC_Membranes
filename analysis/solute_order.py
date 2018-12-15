#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
from LLC_Membranes.llclib import topology, physical
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

        print('Loading trajectory...', end='', flush=True)
        self.t = md.load(traj, top=gro)
        print('Done!')

        self.solute = topology.Solute(solute)
        self.monomer = topology.LC(build_monomer)

        # monomer head group only for testing purposes
        # ref = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
        # solutes_indices = [a.index for a in self.t.topology.atoms if a.residue.name == solute and a.name in ref]
        # self.nsolute = 400
        # self.solute_vectors = self.direction_vectors()
        # self.com = physical.center_of_mass(self.t.xyz[:, solutes_indices, :], [v for v in self.solute.mass.values()][:6])

        solutes_indices = [a.index for a in self.t.topology.atoms if a.residue.name == solute]
        self.nsolute = len(solutes_indices) // self.solute.natoms
        self.solute_vectors = self.direction_vectors()
        self.com = physical.center_of_mass(self.t.xyz[:, solutes_indices, :], [v for v in self.solute.mass.values()])

        reference_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms and
                           a.residue.name in self.monomer.residues]

        self.pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, reference_atoms, :], self.t.unitcell_vectors)

        self.pore_vectors = self.pore_center_vectors()

        self.costheta = self.angles()

        self.nematic_order_parameter = self.nematic()

    def direction_vectors(self):

        print('Calculating direction vector of %s using %s and %s' % (self.solute.name,
                                                                      self.solute.direction_vector[0][0],
                                                                      self.solute.direction_vector[1][0]))

        v = np.zeros([self.t.n_frames, self.nsolute, 2])

        back_atom = [a.index for a in self.t.topology.atoms if a.residue.name == self.solute.name and
                      a.name in self.solute.direction_vector[0]]

        front_atom = [a.index for a in self.t.topology.atoms if a.residue.name == self.solute.name and
                      a.name in self.solute.direction_vector[1]]

        for t in range(self.t.n_frames):
            v[t, ...] = self.t.xyz[t, front_atom, :2] - self.t.xyz[t, back_atom, :2]  # xy projection only

        return v

    def pore_center_vectors(self):

        v = np.zeros([self.t.n_frames, self.nsolute, 2])

        solute_per_pore = self.nsolute // self.npores

        for t in range(self.t.n_frames):
            for p in range(self.npores):
                v[t, (p*solute_per_pore):(p+1)*solute_per_pore, :] = self.pore_centers[t, p, :] - \
                                                             self.com[t, (p*solute_per_pore):(p+1)*solute_per_pore, :2]

        return v

    def angles(self):
        """
        Calculate angle between a vector and plane
        :param vector: vector(s) whose angle with 'plane' we want
        :param plane: vector normal to plane
        :return:
        """

        nT = self.solute_vectors.shape[0]
        costheta = np.zeros([nT, self.nsolute])

        for i in range(nT):
            mag_vsolute = np.linalg.norm(self.solute_vectors[i, ...], axis=1)  # magnitude of each solute vector
            mag_vpore = np.linalg.norm(self.pore_vectors[i, ...], axis=1)
            for s in range(self.nsolute):  # since there is no vectorwise dot product apparently
                costheta[i, s] = np.dot(self.solute_vectors[i, s, :], self.pore_vectors[i, s, :]) / \
                                 (mag_vsolute[s] * mag_vpore[s])

        return costheta

    def nematic(self):
        """
        Calculate the nematic order parameter of a liquid crystal system
        See p. 168 of : http://www.dsf.unica.it/~fiore/libricorsoptr/Chaikin%20Lubensky%20-%20Principles%20of%20Condensed%20Matter%20Physics.jb2.pdf
        :param costheta: cosine of the angle made with director vector (get from "angles" function above)
        :return: nematic order parameter at each frame
        """

        S = np.zeros([self.costheta.shape[0]])  # calculate S for each frame
        for i in range(self.costheta.shape[0]):
            #S[i] = np.sum((3 * self.costheta[i, :-2] ** 2 - 1) / 2) # 3D
            S[i] = np.sum(2 * self.costheta[i, ] ** 2 - 1)  # 2D (https://physics.stackexchange.com/questions/65358/2-d-orientational-order-parameter)

        S /= self.costheta.shape[1]

        return S

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


if __name__ == "__main__":

    args = initialize().parse_args()

    sys = System(args.solute, args.gro, args.traj, args.build_monomer)

    sys.plot(ma=False)
