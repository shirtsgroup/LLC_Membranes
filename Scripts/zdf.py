#! /usr/bin/env python

import numpy as np
import argparse
import mdtraj as md
import tqdm
import matplotlib.pyplot as plt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate a correlation function of positions of atoms in the z direction')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Name of trajectory file (.xtc or .trr)')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop calculations')
    parser.add_argument('--skip', default=1, type=int, help='Usage: --skip n . Sample every nth frame')
    parser.add_argument('-a', '--atoms', nargs='+', type=str, default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'],
                        help='Name of atoms to calculate z distribution function with respect to')

    args = parser.parse_args()

    return args


def z_periodic(positions, unitcell_vectors, images=1):
    """
    :param positions: xyz positions of atoms in original box
    :param unitcell_vectors: the unitcell vectors at each frame so we can calculate z height
    :param images: number of images in the +/- z direction
    :return: new positions with n images in the +/- z direction where n=images
    """

    nT = positions.shape[0]
    natoms = positions.shape[1]
    nimages = images*2 + 1
    periodic = np.zeros([nT, natoms*nimages, 3])
    periodic[:, :natoms, :] = positions
    npores = 4

    # for t in range(nT):
    #     z = np.array([0, 0, box[t, 2, 2]])
    #     # +z direction
    #     for j in range(images):
    #         atom_no = (j + 1)*natoms
    #         for i in range(natoms):
    #             periodic[t, atom_no + i, :] = positions[t, i, :] + (j + 1)*z
    #     for k in range(images):
    #         atom_no = (images + 1 + k)*natoms
    #         for l in range(natoms):
    #             periodic[t, atom_no + l, :] = positions[t, l, :] - (k + 1)*z
    for t in range(nT):
        z = np.array([0, 0, box[t, 2, 2]])
        for p in range(npores):
            for j in range(nimages):


    return periodic


def zdf(positions, npores, atoms_per_layer):
    """
    :param positions: xyz positions of atoms
    :param npores: number of pores
    :param atoms_per_layer: number of atoms in each layer
    :return: number of atoms at a distance z apart
    """

    nT = positions.shape[0]

    natoms = positions.shape[1]
    z = np.zeros([nT, natoms, natoms])  # distance b/w each atom and all other atoms at each frame
    atoms_ppore = natoms / npores
    print natoms, atoms_per_layer, atoms_ppore

    print 'Calculate z distribution function of atom %s' % atom
    for t in tqdm.tqdm(range(nT)):
        for i in range(natoms):
            for j in range(natoms):
                # atoms in the same layer will be at about the same z height. We aren't interested in those values
                # We also only want to calculate distances within pores.
                if int(i/atoms_per_layer) != int(j/atoms_per_layer) and int(i / atoms_ppore) == int(j / atoms_ppore):
                    # ^ i.e. if the atoms aren't in the same layer but they are in the same pore
                    d = positions[t, i, 2] - positions[t, j, 2]
                    z[t, i, j] = d

    z = np.absolute(z).flatten()  # make sure everything is positive and turn it into a 1 dimensional array
    z = np.ma.masked_less(z, 0.1)  # mask the array so all 0 values are not considered. Zero values exist when the z dist
    # is calculated between an atom and itself (i.e. i = j above) or if atoms are in the same layer or different pores

    # bin everything
    bins, edges = np.histogram(z.compressed(), bins=1000)
    bins = [bins[i] / (0.5*nT*natoms*(edges[i + 1] + edges[i])) for i in range(len(bins))]
    # bins = [bins[i] / nT*natoms for i in range(len(bins))]
    plt.figure()
    plt.plot(edges[:-1], bins)
    plt.show()
    exit()

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[args.begin:args.end:args.skip]
    box = t.unitcell_vectors

    for atom in args.atoms:

        keep = [a.index for a in t.topology.atoms if a.name == atom]

        pos = t.atom_slice(keep).xyz

        periodic = z_periodic(pos, box)

        zdf(periodic, 4, 5)