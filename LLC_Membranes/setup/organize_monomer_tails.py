#!/usr/bin/env python

import numpy as np
import mdtraj as md
from LLC_Membranes.llclib import file_rw, topology, transform
import yaml


class Monomer(topology.LC):

    def __init__(self, name, chains, dihedrals):

        super().__init__(name)

        self.nchains = len(chains.keys())
        self.chains = chains
        self.dihedrals = dihedrals
        self.atom_indices = {name: index for index, name in enumerate(self.LC_names)}

        self.organize_bonds()  # makes a dict with keys-value pairs that are atom - atoms bonded to it

    def organize_chains(self):

        for c in range(self.nchains):
            chain = self.chains[c + 1]
            angles = self.dihedrals[c + 1]
            ndihedrals = len(chain) - 3
            for d in range(ndihedrals):
                dih = [[self.atom_indices[i] for i in chain[d:d + 4]]]
                rotation_atoms = self._atom_indices(c + 1, d)
                self.rotate(dih, angles[d], rotation_atoms)

    def _atom_indices(self, chain_no, dihedral_no):
        """ Get all atom indices that need to be rotated along with a given dihedral
        """

        full_chain = [self.atom_indices[a] for a in self.chains[chain_no]]
        atoms = full_chain[(dihedral_no + 3):]

        # This is not general for branched chains. Assumes linear chain with substituents
        rotation_atoms = list(atoms)
        for a in atoms:
            bonded = self.organized_bonds[a]
            for b in bonded:
                if b not in full_chain and b not in rotation_atoms:  # avoid rotated carbon bonded to last dihedral atom
                    rotation_atoms.append(b)

        return rotation_atoms

    def compute_dihedrals(self, dihedral_indices):

        return md.compute_dihedrals(self.t, dihedral_indices, periodic=False) * (180 / np.pi)

    def rotate(self, dihedral_indices, desired_dihedral_angle, rotated_atom_indices):
        """ Determine how many degrees to rotate dihedral

        try following this:
        https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/transformations/rotate.py
        :return:
        """

        current_dihedral = self.compute_dihedrals(dihedral_indices)
        # same as:
        # v1 = np.cross(vectors[0, :], vectors[1, :])
        # v2 = np.cross(vectors[1, :], vectors[2, :])
        # np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))) * (180 / np.pi)
        theta = desired_dihedral_angle - current_dihedral  # angle by which to rotate to achieve desired dihedral

        vectors = np.diff(self.t.xyz[0, dihedral_indices[0], :], axis=0)

        n = vectors[1, :]  # rotation axis

        vector_origin = self.t.xyz[0, dihedral_indices[0][2], :]
        pos = self.t.xyz[0, rotated_atom_indices, :]
        pos -= vector_origin

        R = transform.rotate_about_axis(n, theta)
        rotated = np.zeros_like(pos)
        for i in range(pos.shape[0]):
            rotated[i, :] = np.dot(pos[i, :], R)
        rotated += vector_origin

        self.t.xyz[0, rotated_atom_indices, :] = rotated

    def write_gro(self, out='organized.gro'):

        file_rw.write_gro_pos(self.t.xyz[0, ...], out, ids=self.LC_names, res=self.LC_residues)


if __name__ == "__main__":

    with open('chains.yaml', 'r') as yml:
        cfg = yaml.load(yml)

    lc = Monomer('BOC', cfg['chains'], cfg['angles'])
    lc.organize_chains()
    lc.write_gro()
