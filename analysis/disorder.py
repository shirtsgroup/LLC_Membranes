#!/usr/bin/env python

import argparse
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw, transform
from LLC_Membranes.analysis import Structure_char, Atom_props
from scipy.optimize import minimize
import tqdm
import pickle


def initialize():

    parser = argparse.ArgumentParser(description='Crosslink LLC structure')  # allow input from user

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='wiggle.trr', type=str, help='Name of trajectory to analyze. Preprocess'
                                                                             'so molecules are whole')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Start frame')
    parser.add_argument('-r', '--ref_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Name of '
                        'atoms that will be used to define head groups.')
    parser.add_argument('-pores', default=4, type=int, help='Number of pores in unit cell')
    parser.add_argument('-layers', default=20, type=int, help='Number of monomers per column')
    parser.add_argument('-out', default='disorder.png', type=str, help='')
    parser.add_argument('-pd', '--parallel_displaced', action="store_true", help='Specify if initial configuration is'
                                                                                 'parallel displaced')

    return parser


def angular_deviation(pore, p_center, ideal):

    centered_pore = pore - p_center  # move pore so origin is pore center
    angles = np.arctan2(centered_pore[:, 1], centered_pore[:, 0])  # Calculate angle -- 4 quandrant inverse tangent
    angles_ideal = np.arctan2(ideal[:, 1], ideal[:, 0])

    return angles - angles_ideal


def minimize_deviation(ideal_pore, angle, com, p_center):

    ideal = transform.rotate_coords_z(ideal_pore, angle)  # rotate pore uniformly by some angle

    return np.linalg.norm(angular_deviation(com, p_center, ideal))  # find deviation between ideal and actual positions


class System(object):

    def __init__(self, traj, top, ref_atoms, begin=0, pores=4, layers=20, nrotations=360):

        # initialize what will be the main properties of interest
        self.dz = None
        self.z_values = None
        self.dtheta = None
        self.theta_values = None
        self.dr = None
        self.r_values = None
        self.pore_radius = 0
        self.nrotations = nrotations

        # load trajectory and calculate center of mass of each reference group
        self.t = md.load(traj, top=top)[begin:]

        keep = [a.index for a in self.t.topology.atoms if a.name in ref_atoms]
        natoms = len(ref_atoms)
        self.npores = pores
        self.nlayers = layers
        self.ncol = len(keep) // (natoms * self.nlayers * self.npores)
        self.angle = 360 / self.ncol

        pos = self.t.xyz[:, keep, :]

        ncom = len(keep) // natoms
        self.com = np.zeros([self.t.n_frames, ncom, 3])

        mass = np.array([Atom_props.mass[i] for i in ref_atoms])

        # calculate center of mass of reference groups
        for f in range(self.t.n_frames):
            for i in range(ncom):
                p = pos[f, i * natoms: (i + 1) * natoms, :]
                w = (p.T * mass).T
                self.com[f, i, :] = np.sum(w, axis=0) / sum(mass)

        self.p_centers = Structure_char.avg_pore_loc(self.npores, self.t.xyz, 0)

    def z_deviation(self):

        # map out ideal positions
        zbox = np.mean(self.t.unitcell_vectors[:, 2, 2])
        dbwl = zbox / self.nlayers
        layer_locations = np.linspace(0, dbwl * self.nlayers - dbwl, self.nlayers)

        z = np.zeros([self.t.n_frames, self.com.shape[1]])  # leave out top and bottom layers which are across pbc's
        for f in range(self.t.n_frames):
            for p in range(self.npores):
                for c in range(self.ncol):
                    start = p * self.ncol * self.nlayers + c * self.nlayers
                    end = p * self.ncol * self.nlayers + (c + 1) * self.nlayers
                    z[f, start:end] = self.com[f, start:end, 2] - layer_locations
                # for l in range(1, self.nlayers - 1):
                #     start = p * 100 + l * 5
                #     end = p * 100 + (l + 1) * 5
                #     z[f, start:end] = self.com[f, start:end, 2] - layer_locations[l]

        z = np.where(z >= 0.5*zbox, z-zbox, z)
        z = np.where(z <= -0.5*zbox, z+zbox, z)
        self.z_values = z #.flatten()[np.nonzero(z.flatten())]
        self.z_values -= self.z_values.mean()  # set mean to zero (needed to compare distributions)
        self.dz = np.std(self.z_values)

    def r_deviation(self):

        r = np.zeros([self.t.n_frames, self.com.shape[1]])
        for f in range(self.t.n_frames):
            for p in range(self.npores):
                r[f, p * self.ncol * self.nlayers:(p + 1) * self.ncol * self.nlayers] = np.linalg.norm(
                    self.com[f, p * self.ncol * self.nlayers: (p + 1) * self.ncol * self.nlayers, :2] -
                                                             self.p_centers[:, p, f], axis=1)

        self.pore_radius = np.mean(r.flatten())
        self.r_values = r #.flatten()
        self.dr = np.std(self.r_values)

    def theta_deviation(self, pd=False):

        ideal_column = np.zeros([self.nlayers, 3])
        starting_position = np.array([self.pore_radius, 0, 0])
        ideal_column[:, :] = starting_position
        if pd:
            ideal_column[1::2, :] = transform.rotate_coords_z(starting_position[np.newaxis, :], 33.9155266) #self.angle / 2)
        ideal_pore = np.zeros([ideal_column.shape[0] * self.ncol, 3])
        ideal_pore[:self.nlayers, :] = ideal_column
        for i in range(1, self.ncol):
            ideal_pore[i * self.nlayers: (i + 1) * self.nlayers, :] = transform.rotate_coords_z(
                ideal_column, i * self.angle)

        # layer based systems
        # ideal_layer = np.zeros([5, 3])
        # ideal_layer[0, :] = [self.pore_radius, 0, 0]
        # for i in range(1, 5):
        #     ideal_layer[i, :] = transform.rotate_coords_z(ideal_layer[np.newaxis, 0, :], i * self.angle)
        #
        # ideal_pore = np.zeros([ideal_layer.shape[0] * self.nlayers, 3])
        # for i in range(self.nlayers):
        #     if pd and i % 2 == 1:
        #         ideal_pore[i * ideal_layer.shape[0]: (i + 1) * ideal_layer.shape[0], :] = transform.rotate_coords_z(
        #                                 ideal_layer, self.angle / 2)
        #     else:
        #         ideal_pore[i * ideal_layer.shape[0]: (i + 1) * ideal_layer.shape[0], :] = ideal_layer

        print('Calculating angles needed to minimize angular deviation')
        #angles = np.zeros([self.t.n_frames, self.npores])
        angles = np.linspace(0, 360, self.nrotations)
        angular_deviations = np.zeros([self.t.n_frames, self.npores, self.ncol * self.nlayers])
        for f in tqdm.tqdm(range(self.t.n_frames)):
            for p in range(self.npores):
                deviations = np.zeros([self.nrotations])
                pos = self.com[f, p * self.ncol * self.nlayers:(p + 1) * self.ncol * self.nlayers, :2]
                center = self.p_centers[:, p, f]
                for i, x in enumerate(angles):
                    deviations[i] = minimize_deviation(ideal_pore, x, pos, center)  # , com, ideal_pore, p_centers)
                ideal = transform.rotate_coords_z(ideal_pore, angles[np.argmin(deviations)])
                angular_deviations[f, p, :] = angular_deviation(pos, center, ideal)

        self.theta_values = angular_deviations
        self.theta_values[np.where(self.theta_values < 3)] += 2 * np.pi
        self.theta_values[np.where(self.theta_values > 3)] -= 2 * np.pi
        self.theta_values -= self.theta_values.mean()

        self.dtheta = np.std(self.theta_values)

    def plot_results(self, bins=50, show=True, save=True, out='disorder.png'):

        fig, ax = plt.subplots(1, 3, figsize=(12, 5))

        ax[0].hist(self.z_values.flatten()[np.nonzero(self.z_values.flatten())], bins=bins)
        ax[0].set_title('$\sigma=%.4f nm$' % self.dz)
        ax[0].set_xlabel('COM z distance from ideal layer position (nm)')
        ax[0].set_ylabel('Frequency')

        ax[1].hist(self.r_values.flatten(), bins=bins)
        ax[1].set_title('$\sigma=%.4f nm$' % self.dr)
        ax[1].set_xlabel('COM r distance from pore center (nm)')
        # ax[1].set_ylabel('Frequency')

        ax[2].hist(self.theta_values.flatten(), bins=bins)
        ax[2].set_title('$\sigma=%.4f nm$' % self.dtheta)
        ax[2].set_xlabel('COM $\Theta$ from ideal position (radians)')
        # ax[2].set_ylabel('Frequency')

        plt.tight_layout()
        if save:
            plt.savefig(out)
        if show:
            plt.show()


if __name__ == "__main__":

    args = initialize().parse_args()

    sys = System(args.traj, args.gro, args.ref_atoms, begin=args.begin, pores=args.pores, layers=args.layers)
    sys.z_deviation()
    sys.r_deviation()
    sys.theta_deviation(pd=args.parallel_displaced)

    print(sys.z_values.mean())
    print(sys.theta_values.mean())
    print(sys.r_values.mean())
    # with open('disorder_offset', 'wb') as f:
    #     pickle.dump([sys], f)

    sys.plot_results()
