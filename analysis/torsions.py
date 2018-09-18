#!/usr/bin/env python

import argparse
import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

"""
Calculate the distribution of torsions for a given type of dihedral. Improper dihedrals not implemented! They shouldn't
be moving anyways.
"""


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='traj_whole.xtc', type=str, help='Trajectory file. Molecules should be'
                                                                                'made whole (gmx trjconv -pbc whole)')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-top', '--topology', default='NAcarb11V.itp', type=str, help='Name of topology file for '
                        'molecule that contains torsions of interest')
    parser.add_argument('-d', '--dihedral', nargs='+', default=['C1', 'C2', 'O2', 'C35'], help='Names of atoms, in'
                        'order of their connectivity, that define the dihedral/torsion we are interested in. NOTE: all'
                        'dihedrals of that type will be analyzed')
    parser.add_argument('-r', '--residue', default='HII', type=str, help='Name of residue that dihedral is apart of, as'
                                                                         'it appears in the .gro file')
    parser.add_argument('-exclude', nargs='+', help='Names of atoms to exclude. If that atom appears in a dihedral,'
                                                    'then that dihedral will not be calculated')

    args = parser.parse_args()

    return args


script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # location of this script


def load_topology(top):

    try:
        f = open('%s' % top, 'r')
    except FileNotFoundError:
        try:
            f = open('%s/../top/topologies/%s' % (script_location, top), 'r')
        except FileNotFoundError:
            print('No topology %s found' % top)
            exit()

    monomer_topology = []
    for line in f:
        monomer_topology.append(line)

    f.close()

    return monomer_topology


def calculate_dihedral(pos, indices):
    """
    Follows this post: https://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
    Implements above procedure using numpy for fast vector subtraction, normalization and cross products
                         atom4
                        /
         atom2 --- atom3
        /
    atom1

    :param pos: xyz coordinates of all atoms np.array([nframes, natoms, 3])
    :param indices: indices of atoms forming dihedral for all residues np.array([nresidues, 4]). The 4 indices should
    be ordered according to their connectivity as illustrated above.
    :return: angles of all dihedrals at all frames
    """

    nframes = pos.shape[0]

    b1 = pos[:, indices[:, 1], :] - pos[:, indices[:, 0], :]  # vector from atom1 to atom2 for all dihedrals + frames
    b2 = pos[:, indices[:, 2], :] - pos[:, indices[:, 1], :]  # vector from atom2 to atom3 for all dihedrals + frames
    b3 = pos[:, indices[:, 3], :] - pos[:, indices[:, 2], :]  # vector from atom3 to atom4 for all dihedrals + frames

    # find vectors perpendicular the planes defined by b1 + b2 and b2 + b3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    # normalize
    for t in range(nframes):
        n1[t, :, :] /= np.linalg.norm(n1[t, :, :], axis=1)[:, None]
        n2[t, :, :] /= np.linalg.norm(n2[t, :, :], axis=1)[:, None]
        b2[t, :, :] /= np.linalg.norm(b2[t, :, :], axis=1)[:, None]

    m1 = np.cross(n1, b2)

    # dot product of n1 and n2. faster than np.dot since you have to do it one vector at a time that way
    x = np.zeros([nframes, indices.shape[0]])
    x += n1[:, :, 0] * n2[:, :, 0]
    x += n1[:, :, 1] * n2[:, :, 1]
    x += n1[:, :, 2] * n2[:, :, 2]

    # dot product of m1 and n2
    y = np.zeros([nframes, indices.shape[0]])
    y += m1[:, :, 0] * n2[:, :, 0]
    y += m1[:, :, 1] * n2[:, :, 1]
    y += m1[:, :, 2] * n2[:, :, 2]

    return np.arctan2(y, x) * (180 / np.pi)  # convert to degrees


def ryckaert_belleman(theta, C, mode='degree', normalize=False):
    """
    Evaluate the potential energy of a dihedral at given angle for a potential defined by ryckaert-belleman parameters
    :param theta: angle at which to evaluate potential energy
    :param C: ryckaert_belleman parameters, list of form [C0, C1, C2, C3, C4, C5]
    :param mode: units of theta that are passed to this function. The only option that will do anything is 'degree'
                 since we want everything in radians
    :param normalize: normalize so that max is 1. Useful if plotting on top of angle distribution histogram
    :return: potential energy, float
    """

    theta = np.copy(theta)  # need to make a copy or it will modify input array globally

    if mode == 'degree':
        theta *= (np.pi / 180)

    V = np.zeros_like(theta)
    for i in range(theta.shape[0]):
        for j in range(6):
            V[i] += C[j] * np.cos(theta[i] - np.pi) ** j

    if normalize:
        V /= np.amax(V)

    return V


def largest_prime_factor(n):
    i = 2
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
    return n


def estimated_autocorrelation(x):

    n = len(x)
    variance = x.var()
    # x = x - x.mean()
    r = np.correlate(x, x, mode='full')[-n:]
    result = r / (variance * (np.arange(n, 0, -1)))
    # result = r / np.arange(n, 0, -1)
    return result / np.amax(result[:10])


def autocorrelation_slow(d):
    """
    :param d: numpy array with dihedral angle vs. time
    :return: autocorrelation
    """

    # Subtract mean
    # d -= d.mean(axis=0)

    autocorr = np.zeros([len(d)])
    for l in range(d.shape[0]):
        N = d.shape[0] - l
        for n in range(N):
            autocorr[l] += d[n] * d[n + l]
        autocorr[l] /= N

    autocorr /= d.var()

    return autocorr


class Dihedral(object):

    def __init__(self, gro, traj, top, d, resname, exclusions=[]):
        """
        :param gro: coordinate file (GROMACS .gro format)
        :param traj: trajectory file (GROMACS .trr or .xtc)
        :param top: name of topology file (GROMACS top format). File should be stored in directory where this script
        is located or with all other monomer topologies.
        :param d: names of atoms, in order of connectivity, that define the dihedral/torsion of interest
        :param exclusion: list of atoms to be excluded in dihedral calculations. If that atom is present in exclusions,
        any dihedrals that is a part of will not be counted.
        """

        print('Loading trajectory...', end='', flush=True)
        t = md.load(traj, top=gro)
        print('Done!')

        self.pos = t.xyz
        self.time = t.time
        self.autocorr_fxn = None

        monomer_topology = load_topology(top)

        atoms_index = 0
        while monomer_topology[atoms_index].count('[ atoms ]') == 0:
            atoms_index += 1

        atoms_index += 2  # skip over directive and line of comments (this can be made smarter if necessary)

        types = {}  # atom names with corresponding type
        names = {}  # atom numbers (starting at one) with corresponding names
        natoms_residue = 0  # number of atoms in the residue
        while monomer_topology[atoms_index] != '\n':
            linedata = monomer_topology[atoms_index].split()
            types[linedata[4]] = linedata[1]
            names[linedata[0]] = linedata[4]
            atoms_index += 1
            natoms_residue += 1

        self.dihedral_type = [types[d[0]], types[d[1]], types[d[2]], types[d[3]]]  # we want all dihedrals of this type
        print('Dihedral Type: %s--%s--%s--%s' % (self.dihedral_type[0], self.dihedral_type[1], self.dihedral_type[2],
                                                 self.dihedral_type[3]))

        excluded_atom_numbers = []
        if exclusions:
            for number, name in names.items():
                if name in exclusions:
                    excluded_atom_numbers.append(number)

        # read all the dihedrals
        dihedrals_index = 0
        while monomer_topology[dihedrals_index].count('[ dihedrals ]') == 0:
            dihedrals_index += 1

        dihedrals_index += 2  # skip over directive and line of comments (this can be made smarter if necessary)

        dihedrals = []
        while monomer_topology[dihedrals_index] != '\n':
            dihedral = monomer_topology[dihedrals_index].split()[:4]
            dtype = [types[names[dihedral[i]]] for i in range(4)]
            if dtype == self.dihedral_type or dtype[::-1] == self.dihedral_type:
                # print(dihedral)
                # print(set(dihedral).isdisjoint(exclusions))
                if set(dihedral).isdisjoint(excluded_atom_numbers):  # check that the dihedrals do not contain an excluded atom
                    dihedral = np.array([int(x) for x in dihedral])
                    dihedrals.append(dihedral - 1)  # convert from serial to indices by subtracting 1
            dihedrals_index += 1

        self.ndihedrals = len(dihedrals)  # number of dihedrals to analyze

        res = [a.residue.name for a in t.topology.atoms]
        nres = res.count(resname) // natoms_residue

        self.all_dihedral_angles = np.zeros([self.ndihedrals, t.n_frames, nres])  # weird way to set this up but oh well

        for i in range(self.ndihedrals):
            indices = np.zeros([nres, 4], dtype=int)
            for j in range(nres):
                indices[j, :] = dihedrals[i]*(j + 1)
            self.all_dihedral_angles[i, ...] = calculate_dihedral(self.pos, indices)

    def plot_histogram(self, bins=100, rb=False, ff='gaff', save=False, out='Dihedral.png', show=False):

        hist, bin_edges = np.histogram(self.all_dihedral_angles, bins=bins)
        bin_width = bin_edges[1] - bin_edges[0]
        bin_centers = np.array([bin_edges[i] + bin_width / 2 for i in range(len(bin_edges) - 1)])
        hist = hist / np.amax(hist)  # hist /= np.amax(hist) doesn't work

        if rb:

            forcefield = []
            with open('%s/../top/Forcefields/%s/ffbonded.itp' % (script_location, ff), 'r') as f:
                for line in f:
                    forcefield.append(line)

            dihedral_parameter_index = 0
            while forcefield[dihedral_parameter_index].count('[ dihedraltypes ]') == 0:
                dihedral_parameter_index += 1

            dihedral_parameter_index += 2
            d = forcefield[dihedral_parameter_index].split()[:4]
            while d != self.dihedral_type and d[::-1] != self.dihedral_type:
                dihedral_parameter_index += 1
                d = forcefield[dihedral_parameter_index].split()[:4]
                if dihedral_parameter_index == len(forcefield) - 1:
                    print('Parameter not found >:(')
                    exit()

            C = [float(x) for x in forcefield[dihedral_parameter_index].split()[5:]]

            plt.plot(bin_centers, ryckaert_belleman(bin_centers, C, mode='degree', normalize=True), '--', color='black',
                     label='Dihedral Potential')

        plt.bar(bin_centers, hist, bin_width)
        plt.xlabel('Angle ($\degree$)', fontsize=14)
        plt.ylabel('Count', fontsize=14)
        plt.title('Dihedral Type: %s--%s--%s--%s' % (self.dihedral_type[0], self.dihedral_type[1], self.dihedral_type[2],
                                                     self.dihedral_type[3]))
        plt.tight_layout()
        if save:
            plt.savefig('%s.png' % out)
        if show:
            plt.show()

    def autocorrelation(self, largest_prime=500, cos=False):
        """ FFT based autocorrelation function, which is faster than numpy.correlate. Efficiency is key in order to avoid
        headaches.
        :param largest_prime : the largest prime factor of array length allowed. The smaller the faster. 1.6M points takes
        about 5 seconds with largest_prime=1000. Just be aware that you are losing data by truncating. But 5-6 data points
        isn't a big deal for large arrays.
        """

        shape = self.all_dihedral_angles.shape
        # np.reshape won't quite work how I need it to. This way is ~50x slower than reshape but still fast enough

        dihedrals = np.zeros([shape[1], shape[0]*shape[2]])  # of size [nframes, ndihedrals]

        for i in range(shape[0]):
            dihedrals[:, i*shape[2]:(i+1)*shape[2]] = self.all_dihedral_angles[i, ...]

        if cos:
            dihedrals *= (np.pi/180)
            dihedrals = np.cos(dihedrals)

        # Don't allow a prime factor larger than 'largest_prime'. Truncate data until that condition is met
        l = 2 * dihedrals.shape[0] - 1

        while largest_prime_factor(l) >= largest_prime or l % 2 == 0:
            l -= 1

        dihedrals = dihedrals[:(l + 1) // 2, :]
        length = dihedrals.shape[0] * 2 - 1

        dihedrals -= np.mean(dihedrals, axis=0)

        fftx = np.fft.fft(dihedrals, n=length, axis=0)
        ret = np.fft.ifft(fftx * np.conjugate(fftx), axis=0)
        ret = np.fft.fftshift(ret, axes=(0,))

        self.autocorr_fxn = ret[length // 2:].real
        self.autocorr_fxn /= np.arange(dihedrals.shape[0], 0, -1)[:, None]
        self.autocorr_fxn /= np.var(dihedrals, axis=0)

        self.autocorr_fxn = np.mean(self.autocorr_fxn, axis=1)

        return self.autocorr_fxn  # normalized

    def plot_autocorrelation(self, show=False, save=True, savename='dihedral_autocorrelation'):

        plt.figure()
        plt.plot(self.time/1000, self.autocorr_fxn, linewidth=2)
        #plt.plot([self.time[0]/1000, self.time[-1]/1000], [0, 0], '--', color='black')
        plt.xlabel('Time (ns)', fontsize=14)
        plt.ylabel('Autocorrelation', fontsize=14)
        plt.tight_layout()
        if save:
            plt.savefig('%s.pdf' % savename)
        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize()

    dihedrals = Dihedral(args.gro, args.traj, args.topology, args.dihedral, args.residue, exclusions=args.exclude)
    dihedrals.plot_histogram(rb=True, save=True, show=True)
    dihedrals.autocorrelation(cos=True)
    dihedrals.plot_autocorrelation(show=True)
