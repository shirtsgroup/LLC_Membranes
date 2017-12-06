#! /usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from llclib import physical
import place_solutes


def initialize():

    parser = argparse.ArgumentParser(description='Run Cylindricity script')

    parser.add_argument('-t', '--traj', default='wiggle.trr', help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file')
    parser.add_argument('-b', '--begin', default=0, type=int, help='End frame')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Start frame')
    parser.add_argument('-ref', default='C6', type=str, help='Reference atom for locating pore centers')
    parser.add_argument('-cut', default=1, type=float, help='Distance from pore water molecules where water molecules'
                                                            'will no be removed')
    parser.add_argument('-gap', action="store_true", help='Use this flag if there is a gap. Water in the gap will not'
                        'be considered while searching for hydrogen bonds. Do not use sodium as the reference atom in '
                        'this case since ions will diffuse up into the water')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    print('Loading trajectory...', end="")
    t = md.load('%s' % args.traj, top='%s' % args.gro)[args.begin:args.end]
    print('Done!')

    npores = 4

    ref = [a.index for a in t.topology.atoms if a.name == args.ref]

    p_centers = physical.avg_pore_loc(npores, t.xyz[:, ref, :])

    # get indices of all water molecule 'O' atoms. Restrict to 'O' so parts of water molecules won't get cut off later
    w = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']

    water = t.xyz[:, w, :]

    blacklist = []  # all the water molecules which will be removed

    # if args.gap:
    #     # remove all water above and below the membrane using the last frame as a reference
    #     top = np.max(t.xyz[-1, ref, 2])  # max of z positions of all reference atoms in last frame
    #     bot = np.min(t.xyz[-1, ref, 2])  # min of z positions of all reference atoms in last frame
    #     for i in range(water.shape[1]):
    #         if top < water[-1, i, 2] or water[-1, i, 2] < bot:
    #             blacklist.append(w[i])

    # remove all water molecules within args.cut of the pore centers. Use just the last frame to make this decision
    for p in range(npores):
        # narrow down the positions to those that are with 'cut' of at least one pore
        distances = np.linalg.norm(water[-1, :, :2] - p_centers[:, p, -1], axis=1)

        for i, d in enumerate(distances):
            if d < args.cut:  # if the water molecule is located too close to the pore center
                if w[i] not in blacklist:
                    blacklist.append(w[i])  # add OW to the blacklist for removal
                    blacklist.append(w[i] + 1)  # HW1 attached to OW
                    blacklist.append(w[i] + 2)  # HW2 attached to OW

    keep = [a.index for a in t.topology.atoms if a.index not in blacklist]

    print("Getting rid of unwanted water molecules...", end="")
    mod = t.atom_slice(keep)
    print("Done!")

    print("Calculating hbonds...", end="")
    hbonds = md.baker_hubbard(mod, periodic=True, exclude_water=False)
    print("Done!")

    print(hbonds)
    # label = lambda hbond : '%s -- %s' % (t.topology.atom(hbond[0]), t.topology.atom(hbond[2]))
    #
    # for hbond in hbonds:
    #     print(label(hbond))

    # from llclib import file_rw
    # full_box = t.unitcell_vectors
    #
    # box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
    #                full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]
    #
    # ids = [a.name for a in t.topology.atoms if a.index not in blacklist]
    # res_names = [a.residue.name for a in t.topology.atoms if a.index not in blacklist]
    # file_rw.write_gro_pos(t.xyz[-1, keep, :], 'rm.gro', box=box_gromacs, ids=ids, res=res_names)
    #
    #     hist, bin_edges = np.histogram(d_sorted[:stop], bins=nbins, range=(0, cut))  # the range option is necessary
    #     #  to make sure we have equal sized bins on every iteration
    #
    #     density += hist / box[t, 2, 2]
    #
