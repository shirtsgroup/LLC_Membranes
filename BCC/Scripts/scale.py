#!/usr/bin/env python

import mdtraj as md
import numpy as np
import argparse
import bcc_class
from llclib import transform
from llclib import file_rw


def initialize():

    parser = argparse.ArgumentParser(description='Isotropically scale a bcc structure')  # allow input from user

    parser.add_argument('-g', '--gro', default='Ordered.gro', type=str, help='Name of input file')
    parser.add_argument('-o', '--output', default='scaled.gro', help='Name of reordered output file')
    parser.add_argument('-m', '--monomer', default='Dibrpyr14.gro', help='Name of monomer used to construct structure')
    parser.add_argument('-f', '--factor', type=float, default=2.0, help='Factor to scale coordinates by')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)  # load all information from .gro file
    pos = t.xyz[0, :, :]

    LC = bcc_class.LC(args.monomer)  # create a liquid crystal object for the chosen monomer

    ref = [a.index for a in t.topology.atoms if a.name in LC.ref_names]  # get the indices of the reference atoms defined by the class
    nref = len(LC.ref_names)  # number of reference atoms

    ref_slice = t.atom_slice(ref).xyz  # get positions of all reference atoms
    pos_ref = np.zeros([int(ref_slice.shape[1] / nref), 3])  # we want a single point for each molecule

    for i in range(pos_ref.shape[0]):  # average the positions of the reference atoms within each monomer
        for j in range(nref):
            pos_ref[i, :] += ref_slice[0, i*nref + j, :]

    pos_ref /= (nref / args.factor)  # average coordinates and apply scaling factor

    # Now translate all monomers to the new reference positions
    nres = int(t.n_atoms / LC.natoms)  # number of residues
    non_ion = LC.natoms - LC.no_ions  # number of atoms in residue not including ions. (Ions are grouped as a separate residues)

    scaled = np.zeros([LC.natoms*nres, 3])

    for i in range(nres):
        scaled[i*non_ion:(i + 1)*non_ion, :] = transform.translate(pos[i*non_ion:(i + 1)*non_ion, :], pos_ref[i, :] / args.factor, pos_ref[i, :])

    for i in range(nres):
        for j in range(LC.no_ions):
            scaled[(non_ion*nres + i*LC.no_ions):(non_ion*nres + (i + 1)*LC.no_ions), :] = \
            transform.translate(pos[(non_ion*nres + i*LC.no_ions):(non_ion*nres + (i + 1)*LC.no_ions), :], pos_ref[i, :] / args.factor, pos_ref[i, :])

    res = LC.resid[:non_ion]*nres
    res += LC.resid[non_ion:]*nres
    ids = LC.names[:non_ion]*nres
    ids += LC.names[non_ion:]*nres

    file_rw.write_gro_pos(scaled, '%s' % args.output, res=res, ids=ids, box=args.factor*t.unitcell_lengths[0])
