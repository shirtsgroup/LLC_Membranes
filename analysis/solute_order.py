#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
from LLC_Membranes.llclib import topology, physical


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


class System(object):

    def __init__(self, solute, gro, traj, build_monomer, spline=False, npores=4):

        self.npores = npores

        print('Loading trajectory...', end='', flush=True)
        self.t = md.load(traj, top=gro)
        print('Done!')

        self.solute = topology.Solute(solute)
        self.monomer = topology.LC(build_monomer)

        solutes_indices = [a.index for a in self.t.topology.atoms if a.residue.name == solute]
        self.nsolute = len(solutes_indices) // self.solute.natoms
        self.solute_vectors = self.direction_vectors()
        self.com = physical.center_of_mass(self.t.xyz[:, solutes_indices, :], [v for v in self.solute.mass.values()])

        reference_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms and
                           a.residue.name in self.monomer.residues]

        self.pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, reference_atoms, :])

        self.pore_vectors = self.pore_center_vectors()
        print(self.pore_vectors.shape)

    def direction_vectors(self):

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

        for t in range(self.t.n_frames):
            for p in range(self.npores):
                v[t, (p*self.npores):(p+1)*self.npores, :] = self.com[t, (p*self.npores):(p+1)*self.npores, :2] - \
                                                            self.pore_centers[t, p, :]

        return v


if __name__ == "__main__":

    args = initialize().parse_args()

    sys = System(args.solute, args.gro, args.traj, args.build_monomer)
