#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
from builtins import range
from past.utils import old_div
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
    parser.add_argument('--avg', action="store_true", help='Average the distributions of all selected atoms')
    parser.add_argument('-bins', default=1000, type=int, help='Number of bins to use when binning distances')
    parser.add_argument('-apl', default=5, type=float, help='Atoms or monomers per layer')

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
    npores = 4
    atomppore = old_div(natoms, npores)

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
    m = np.linspace(-images, images, nimages)
    for t in range(nT):
        z = np.array([0, 0, box[t, 2, 2]])
        for p in range(npores):
            for j in range(nimages):
                for k in range(atomppore):
                    periodic[t, nimages*atomppore*p + atomppore*j + k, :] = positions[t, atomppore*p + k, :] + m[j]*z

    return periodic


def zdf(positions, npores, atoms_per_layer, box, images=1, tol=0.01):
    """
    :param positions: xyz positions of atoms expanded periodically in the +/- z direction
    :param npores: number of pores
    :param atoms_per_layer: number of atoms in each layer
    :return: number of atoms at a distance z apart
    """

    nT = positions.shape[0]
    nimages = images*2 + 1
    natoms = positions.shape[1]

    z = np.zeros([nT, natoms, natoms])  # distance b/w each atom and all other atoms at each frame
    atoms_ppore = old_div(natoms, npores)

    print('Calculating z distribution function of atom %s' % atom)
    for t in tqdm.tqdm(list(range(nT))):
        for p in range(npores):

            # ensure we are only looking at pores from the original unit cell (the middle image in the periodic cell created by z_periodic)
            start = p*atoms_ppore + (old_div(atoms_ppore, nimages))
            end = start + (old_div(atoms_ppore, nimages))
            for i in range(start, end):
                pore_start = p * atoms_ppore
                pore_end = (p + 1) * atoms_ppore
                for j in range(pore_start, pore_end):
                    # now calculate distance between atoms in original unit cell and all atoms in periodic cell
                    # atoms in the same layer will be at about the same z height. We aren't interested in those values
                    # We also only want to calculate distances within pores.
                    if int(old_div(i,atoms_per_layer)) != int(old_div(j,atoms_per_layer)): # and int(i / atoms_ppore) == int(j / atoms_ppore):
                        # ^ if the atoms aren't in the same layer but they are in the same pore ^
                        d = positions[t, i, 2] - positions[t, j, 2]
                        z[t, i, j] = d

    L = np.mean(box[:, 2, 2])  # average z length of unit cell

    p = old_div(atoms_ppore, ((images*2 + 1)*L))  # average number density of atoms [particles / nm]

    # in periodic cell. Atoms in same layer are not counted (particles / nm)

    z = np.absolute(z).flatten()  # make sure everything is positive and turn it into a 1 dimensional array
    z = np.ma.masked_equal(z, 0)  # mask the array so all 0 values are not considered. Zero values exist when the z dist
    # is calculated between an atom and itself (i.e. i = j above) or if atoms are in the same layer or different pores

    # bin everything
    bins, edges = np.histogram(z.compressed(), bins=args.bins)

    end = 0
    while edges[end] < L:
        end += 1

    dx = edges[-1] - edges[-2]

    bins = [old_div(i, (2 * nT * p * dx * L * npores)) for i in bins]  # normalize by number density

    # control = np.zeros([len(bins)]) + np.mean(bins[:end])

    end /= 2

    if not args.avg:
        plot_zdf(edges, bins, end=end)

    return edges[:int(end)], bins[:int(end)]


def plot_zdf(z, d, end=-1):
    """
    :param z: distance along membrane z dimension ([nbins])
    :param d: density at each z location
    :return: plot the zdf
    """

    plt.figure()
    plt.plot(z[:end], d[:end])
    # plt.title('Z distribution function')
    plt.xlabel('Z distance separation (nm)', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    plt.savefig('zdf.png')
    plt.show()


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[args.begin:args.end:args.skip]
    frame = t.slice(0)
    # from llclib import file_rw
    # file_rw.write_gro(frame, 'first.gro')
    # exit()
    box = t.unitcell_vectors
    L = np.mean(box[:, 2, 2])  # average z length of unit cell

    if args.avg:

        zdf_avg = np.zeros([args.bins])

    for atom in args.atoms:

        keep = [a.index for a in t.topology.atoms if a.name == atom]

        pos = t.atom_slice(keep).xyz

        periodic = z_periodic(pos, box)

        # from llclib import file_rw
        # file_rw.write_gro_pos(periodic[0, :, :], 'first')
        # exit()

        z, d = zdf(periodic, 4, args.apl, box)

        if args.avg:

            zdf_avg[:len(d)] += d

    if args.avg:

        zdf_avg /= len(args.atoms)  # average of all zdfs

        zdf_avg = np.trim_zeros(zdf_avg, trim='b')

        plot_zdf(z, zdf_avg[:z.shape[0]])