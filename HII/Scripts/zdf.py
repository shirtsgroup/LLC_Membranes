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
    parser.add_argument('-l', '--load', type=str, help='Load compressed numpy array')
    parser.add_argument('-com', '--com', action="store_true", help='Calculate zdf based on com of atoms')
    parser.add_argument('-r', '--res', type=str, help='Residue name to calculate zdf of')

    args = parser.parse_args()

    return args


def centroids(x, n):
    """
    Calculate the centroid of group of atoms based on atomic positions.
    Same thing as center of mass if all atoms are the same
    :param x: xyz positions of all atoms
    :param n: number of atoms in each group (assumed atoms are grouped sequentially)
    :return: centroids
    """

    nT = x.shape[0]
    natoms = x.shape[1]
    ngrps = int(natoms / n)

    centers = np.zeros([nT, ngrps, 3])

    for t in range(nT):
        for i in range(ngrps):
            for j in range(int(n)):
                centers[t, i, :] += x[t, i*int(n) + j, :]

    centers /= n

    return centers


def z_periodic(positions, box, images=1):
    """
    :param positions: xyz positions of atoms in original box
    :param box: the unitcell vectors at each frame so we can calculate z height
    :param images: number of images in the +/- z direction
    :return: new positions with n images in the +/- z direction where n=images
    """

    # use cKDTree to make this way faster
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


def zdf(positions, npores, atoms_per_layer, box, bins, name, images=1, tol=0.01):
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

    print('Calculating z distribution function of atom %s' % name)
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
    bins, edges = np.histogram(z.compressed(), bins=bins)

    end = 0
    while edges[end] < L:
        end += 1

    dx = edges[-1] - edges[-2]

    bins = [old_div(i, (2 * nT * p * dx * L * npores)) for i in bins]  # normalize by number density

    # control = np.zeros([len(bins)]) + np.mean(bins[:end])

    end /= 2
    # print(edges, bins)
    # exit()
    # if not args.avg:
    #     plot_zdf(edges, bins, end=end)

    return edges[:int(end)], bins[:int(end)]


def plot_zdf(z, d, out, end=-1):
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
    plt.axes().tick_params(labelsize=14)
    plt.tight_layout()
    plt.savefig('%s' % out)
    plt.show()


def spacing(x, y, start=0.3, window=0.25, zcut=4):
    """
    Find the layer spacing based on the first maximum of the zdf
    :param x: z values where zdf is taken
    :param y: zdf values for each z
    :param dmax: maximum distance between layers. Search for maxima stops at this distance
    :return: layer spacing
    """

    bin_width = x[1] - x[0]
    window_bins = int(window / bin_width)

    begin = 0
    while x[begin] < window:
        begin += 1
    begin -= 1

    cut = 0
    while x[cut] < zcut:
        cut += 1

    maxes = [[], []]
    mins = [[], []]
    max_search = 1
    while (begin + window_bins) < cut:

        if max_search == 1:
            m = np.max(y[begin:(begin+window_bins)])
            max_search = 0
            m_ndx = np.where(y == m)[0][0]
            maxes[1].append(m)
            maxes[0].append(x[m_ndx])
        elif max_search == 0:
            m = np.min(y[begin:(begin+window_bins)])
            max_search = 1
            m_ndx = np.where(y == m)[0][0]
            mins[1].append(m)
            mins[0].append(x[m_ndx])

        begin += (m_ndx - begin)

    separation = []
    fluct = []
    separation.append(maxes[0][0])
    for i in range(1, len(maxes[1])):
        separation.append(maxes[0][i] - maxes[0][i - 1])
        fluct.append(mins[1][i - 1] - maxes[1][i - 1])
    fluct.append(mins[1][-1] - maxes[1][-1])

    return maxes[0][0], abs(fluct[0] / 2)


def power_spectrum(data, bin):
    """
    Compute the power spectrum of the data (find dominant frequencies in the fourier series fitting discrete data)
    :param data: Data to be fourier transformed (1D numpy array)
    :param bin: bin size (nm)
    :return: power spectrum of data (ps) with corresponding frequencies (freqs), the max frequency (max) and
    """

    data = data - np.mean(data)  # get rid of a peak at zero
    ps = np.abs(np.fft.fft(data))**2

    freqs = np.fft.fftfreq(data.size)
    idx = np.argsort(freqs)

    # fft = np.abs(np.fft.fft(data))
    # flat = fft.flatten()
    # flat.sort()
    #
    # max_freq = np.where(fft == flat[-1])[0][0]
    max_freq = np.argmax(np.abs(np.fft.fft(data)))
    freq = freqs[max_freq]
    # if freq == 0:  # sometimes there is a large spike at 0
    #     max_freq = np.where(fft == flat[-2])[0][0]
    #     freq = freqs[max_freq]

    # modify things so they'll plot nicely and in the correct units
    freqs = freqs[idx] / bin
    ps = ps[idx]
    max = abs(freq / bin)  # maximum frequency in hertz

    return ps, freqs, max


if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)[args.begin:args.end:args.skip]
    frame = t.slice(0)

    # from llclib import file_rw
    # file_rw.write_gro(frame, 'first.gro')
    # exit()
    box = t.unitcell_vectors
    L = np.mean(box[:, 2, 2])  # average z length of unit cell

    if args.load:

        zdf = np.load(args.load)
        zdf_avg = zdf["zdf_avg"]
        z = zdf["z"]
        bin_width = (z[1] - z[0])
        ps, freqs, max = power_spectrum(zdf_avg, bin_width)
        # dbwl, amp = spacing(z, zdf_avg)
        print('Fourier distance between layers: %s nm' % (1/max))
        # print('Amplitude of first peak : %2.2f %% of mean' % (100*(amp/np.mean(zdf_avg[:z.shape[0]]))))
        plt.figure()
        positive = 0
        while freqs[positive] < 0:
            positive += 1
        plt.plot(freqs[positive:], ps[positive:] / np.max(ps))
        plt.xlabel('Frequency (cycle/nm)', fontsize=14)
        plt.ylabel('Normalized Intensity', fontsize=14)
        plt.axes().tick_params(labelsize=14)
        plt.tight_layout()
        plt.savefig('ps.png')
        plot_zdf(z[:len(zdf_avg)], zdf_avg[:z.shape[0] - 5])

    else:

        if args.com:

            keep = [a.index for a in t.topology.atoms if a.name in args.atoms]
            pos = t.atom_slice(keep).xyz
            centers = centroids(pos, args.apl)
            periodic = z_periodic(centers, box)
            atom = 'center of mass'
            z, d = zdf(periodic, 4, args.apl, box, args.bins, atom)
            d = np.trim_zeros(d, trim='b')
            plot_zdf(z[:-2], d[:z.shape[0] - 2], 'zdf.png')
            # dbwl = spacing(z, d)
            # print('Distance between layers: %s' % dbwl)
            exit()

        if args.avg:

            zdf_avg = np.zeros([args.bins])

        for atom in args.atoms:

            if args.res:
                keep = [a.index for a in t.topology.atoms if a.name == atom and a.residue.name == args.res]
            else:
                keep = [a.index for a in t.topology.atoms if a.name == atom]

            # pos = t.atom_slice(keep).xyz
            pos = t.xyz[:, keep, :]

            periodic = z_periodic(pos, box)

            # from llclib import file_rw
            # file_rw.write_gro_pos(periodic[0, :, :], 'first')
            # exit()

            z, d = zdf(periodic, 4, args.apl, box, args.bins, atom)

            if args.avg:

                zdf_avg[:len(d)] += d

        if args.avg:

            zdf_avg /= len(args.atoms)  # average of all zdfs

            zdf_avg = np.trim_zeros(zdf_avg, trim='b')

            np.savez_compressed("zdf", zdf_avg=zdf_avg, z=z)

            bin_width = (z[1] - z[0])
            ps, freqs, max = power_spectrum(zdf_avg, bin_width)
            print('Fourier distance between layers: %s nm' % (1/max))
            plt.figure()
            plt.plot(freqs, ps)
            # dbwl = spacing(z, zdf_avg)
            plot_zdf(z, zdf_avg[:z.shape[0]], 'zdf.png')
            # dbwl = spacing(z, zdf_avg)
            # print('Distance between layers: %s' % dbwl)