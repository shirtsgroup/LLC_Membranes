#! /usr/bin/env python

import argparse
import mdtraj as md
from place_solutes import trace_pores
import numpy as np
import matplotlib.pyplot as plt
import tqdm
from scipy import spatial
import Atom_props


def initialize():

    parser = argparse.ArgumentParser(description='Figure out the weight percent of water in the pores and tails')

    parser.add_argument('-t', '--traj', default='PR.xtc', help='Name of GROMACS trajectory file')
    parser.add_argument('-g', '--gro', default='PR.gro', help='Name of GROMACS coordinate file')
    parser.add_argument('-ox', '--tail_oxygen', nargs='+', default=['O5', 'O6', 'O7', 'O8', 'O9', 'O10'], help='Oxygen'
                        'atoms that will be used to define the tail region')
    parser.add_argument('-tr', '--tail_radius', default=0.5, type=float, help='Max distance from tail oxygens a water '
                        'molecule can exist in order to be counted as inside the pore')
    parser.add_argument('-p', '--pore_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms that'
                        'will be used to define the pore region')
    parser.add_argument('-pr', '--pore_radius', default=0.5, type=float, help='Max distance from pore center a water '
                        'molecule can exist in order to be counted as inside the pore')
    parser.add_argument('-b', '--bounds', default=3, type=float, help='Distance from center up until which all atoms '
                                                                      'will be included in calculation (nm)')
    parser.add_argument('-natoms', default=137, type=int, help='Number of atoms in monomer residue (not including ions '
                                                               ' if they are separate residues!')
    parser.add_argument('-mpl', default=5, type=int, help='Number of monomers per layer')
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of layers')
    parser.add_argument('--load', action="store_true")
    parser.add_argument('--save', action="store_true")

    args = parser.parse_args()

    return args


def make_groups(t, pos, pore_atoms, tail_atoms, bounds, natoms, mpl):
    """
    :param t: mdtraj trajectory object
    :param pore_atoms: atoms that define the pore region, list
    :param tail_atoms: atoms that define the tail region, list
    :param bounds: only include atoms within these z boundaries. list:[lower_bound, upper_bound]
    :param natoms: number of atoms in monomer residue
    :param mpl: number of monomers per layer
    :return:
    """

    all_pore_atoms = []
    all_tail_atoms = []
    water = []  # oxygen atom of water molecule
    other = []  # all other atoms within bounds. They will be needed for a weight percent calculation
    mass = 0
    for a in t.topology.atoms:
        if bounds[0] <= pos[a.index, 2] <= bounds[1]:
            mass += Atom_props.mass[a.name]
            if a.residue.name == 'HOH' and a.name == 'O':
                water.append(a.index)
            elif a.name in pore_atoms:
                all_pore_atoms.append(a.index)
            elif a.name in tail_atoms:
                all_tail_atoms.append(a.index)

    # all_pore_atoms will be used to make a spline which traces through the pores. To do that properly, only full layers
    # can be included in the list.

    # layer_number = [i // (mpl*natoms) for i in all_pore_atoms]
    #
    # pore_atoms_keep = []
    # for i in range(len(all_pore_atoms)):
    #     if layer_number.count(layer_number[i]) == int(mpl*len(pore_atoms)):
    #         pore_atoms_keep.append(all_pore_atoms[i])
    #     else:
    #         other.append(all_pore_atoms[i])

    return all_pore_atoms, all_tail_atoms, water, mass


def restrict_spline(full_spline, bounds):
    """
    Restrict a pore-tracing spline to be within two z boundaries
    :param full_spline: A full spline
    :param bounds: z boundaries
    :return: restricted spline to within z-boundaries
    """

    restricted = []
    for i in range(full_spline.shape[0]):
        if bounds[0] <= full_spline[i, 2] <= bounds[1]:
            restricted.append(full_spline[i, :])

    return np.array(restricted)


def water_content(pos, ref_pos, r):
    """
    Find number of water molecules in region
    :param pos: positions of all water molecules (excluding those outside boundaries defined by args.bounds)
    :param ref_pos: positions defining the region
    :param r: distance from reference positions a point can be for it to be included as a part of the region
    :return: number of water molecules in region
    """

    n = 0
    tree = spatial.cKDTree(ref_pos)
    for i in range(pos.shape[0]):
        d = tree.query(pos[i, :])[0]
        if d < r:
            n += 1

    return n


if __name__ == "__main__":

    args = initialize()

    print('Loading trajectory...')
    t = md.load(args.traj, top=args.gro)
    print('Done!')
    pos = t.xyz
    box = t.unitcell_vectors
    res = np.array([a.residue.name for a in t.topology.atoms])
    ids = np.array([a.name for a in t.topology.atoms])
    pores = [a.index for a in t.topology.atoms if a.name in args.pore_atoms]
    box_gromacs = [box[0, 0, 0], box[0, 1, 1], [0, 2, 2], [0, 0, 1], [0, 2, 0],
                   box[0, 1, 0], box[0, 0, 2], [0, 1, 2], [0, 2, 0]]
    mwater = 18.016
    disregard = False

    nframes = pos.shape[0]
    wt_pores = np.zeros([nframes])
    wt_tails = np.zeros([nframes])
    combined = np.zeros([nframes])
    wt_tot = np.zeros([nframes])

    for f in tqdm.tqdm(range(nframes)):

        middle = box[f, 2, 2] / 2  # ~ halfway into the membrane based on last frame

        upper_bound = middle + args.bounds
        lower_bound = middle - args.bounds

        all_pore_atoms, all_tail_atoms, water, mass = make_groups(t, t.xyz[f, :, :], args.pore_atoms, args.tail_oxygen,
                                                                   [lower_bound, upper_bound], args.natoms, args.mpl)

        full_pore_spline = trace_pores(pos[f, pores, :], box[f, :, :], 20)

        spline = restrict_spline(full_pore_spline, [lower_bound, upper_bound])

        n_water_pore = water_content(pos[f, water, :], spline, args.pore_radius)
        n_water_tail = water_content(pos[f, water, :], pos[f, all_tail_atoms, :], args.tail_radius)

        if n_water_pore + n_water_tail > len(water):
            if not disregard:
                print('WARNING: Water in tail and water in pore add up to more than the total amount of water in the '
                      'restricted unit cell. Your pore and tail radii are probably too large.')
                user = input('Would you like to continue and ignore this message (y/n)')
                if user == 'y':
                    disregard = True
                else:
                    exit()

        wt_pores[f] = n_water_pore*mwater / mass
        wt_tails[f] = n_water_tail*mwater / mass
        combined[f] = wt_pores[f] + wt_tails[f]
        wt_tot[f] = len(water)*mwater / mass

    plt.plot(t.time, 100*wt_pores, label='Pores')
    plt.plot(t.time, 100*wt_tails, label='Tails')
    plt.plot(t.time, 100*combined, label='Pores+Tails')
    plt.plot(t.time, 100*wt_tot, label='Total')
    plt.legend()
    plt.xlabel('Time (ps)')
    plt.ylabel('Water content (%)')
    plt.tight_layout()
    plt.savefig('water_content.png')
    plt.show()


