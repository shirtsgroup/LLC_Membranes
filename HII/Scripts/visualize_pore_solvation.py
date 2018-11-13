#!/usr/bin/env python

import numpy as np
import argparse
import mdtraj as md
import subprocess


def initialize():

    parser = argparse.ArgumentParser(description='Remove all atoms except those selected and nearby water molecules')

    parser.add_argument('-t', '--traj', type=str, default='PR.xtc', help='Simulation trajectory')
    parser.add_argument('-g', '--gro', type=str, default='wiggle.gro', help='Coordinate file')
    parser.add_argument('-a', '--atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms that make'
                        'up the pore wall')
    parser.add_argument('-b', '--buffer', type=float, default=0.2, help='Fraction of unitcell to exclude from each side')
    parser.add_argument('-r', '--radius', type=float, default=0.5, help='Max distance a water molecule can be from'
                                                                        'a pore wall atom to not be thrown away')
    parser.add_argument('-o', '--output', type=str, default='frames', help='Output trajectory')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.traj, top=args.gro)
    z = t.unitcell_vectors[-1, 2, 2]  # just use last frame box vector
    buffer = args.buffer * z
    top = z - buffer
    bottom = buffer

    pore_wall_ndx = [a.index for a in t.topology.atoms if a.name in args.atoms]

    all_water_ndx = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' if bottom < t.xyz[-1, a.index, 2] < top]

    pore_wall = t.xyz[-1, pore_wall_ndx, :]

    water_to_keep = []
    for i in range(len(all_water_ndx)):
        distances = np.linalg.norm(pore_wall - t.xyz[-1, all_water_ndx[i], :], axis=1)
        if np.amin(distances) < args.radius:
            water_to_keep.append(all_water_ndx[i])

    atoms_to_keep = pore_wall_ndx + water_to_keep

    with open('pore_solvation.tcl', 'w') as f:

        f.write('mol modselect 0 0 index ')
        for i in pore_wall_ndx:
            f.write('%s ' % i)
        f.write('\n')
        f.write('mol addrep 0\n')
        f.write('mol modselect 1 0 index ')
        for i in water_to_keep:
            f.write('%s ' % i)

    subprocess.Popen(["vmd", args.gro, args.traj, '-e', 'pore_solvation.tcl'])
