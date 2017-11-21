#!/usr/bin/env python

import argparse
import mdtraj as md
import place_solutes
import lc_class
import numpy as np
import top
from llclib import file_rw


def initialize():

    parser = argparse.ArgumentParser(description='Place solutes in the pores')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file where solutes will be placed')
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of layers in input system')
    parser.add_argument('-b', '--buffer', default=2, type=int, help='Number of layers to remove on each side of the '
                        'membrane')
    parser.add_argument('-r', '--ref', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Reference atom (s)'
                        'used to locate pores')
    parser.add_argument('-m', '--build_monomer', default='NAcarb11Vd', help='Name of monomer used to build system')
    parser.add_argument('-o', '--out', type=str, default='stitched.gro', help='Name of output .gro file')
    parser.add_argument('-c', '--configuration', type=str, default='offset', help='Orientation of aromatic head groups')
    parser.add_argument('--charge', type=float, default=1, help='Magnitude of charge to place on vsites')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    t = md.load(args.gro)
    full_box = t.unitcell_vectors
    monomer = lc_class.LC('%s.gro' % args.build_monomer)

    # remove buffered layers
    nmon = -2  # subtract two for the top lines in the .gro
    ids_full = []
    res_full = []
    with open(args.gro) as f:
        for line in f:
            if str.strip(line[5:10]) in monomer.ions:
                break
            else:
                nmon += 1
            if nmon > 0:
                ids_full.append(str.strip(line[10:15]))  # will need these to create .gro later
                res_full.append(str.strip(line[5:10]))

    nmon //= monomer.natoms  # double divide makes the output an integer (floors it too)
    mon_per_layer = nmon // 4 // args.layers
    mon_atoms_per_pore = mon_per_layer*args.layers*monomer.natoms

    # get indices of atoms to keep (i.e. don't include buffer layers)
    mon = []
    ions = []
    ids = []
    res = []
    for p in range(4):
        for l in range(args.buffer, args.layers - args.buffer):
            for m in range(mon_per_layer):
                for a in range(monomer.natoms):
                    ndx = p*mon_atoms_per_pore + l*mon_per_layer*monomer.natoms + m*monomer.natoms + a
                    mon.append(ndx)
                    ids.append(ids_full[ndx])
                    res.append(res_full[ndx])

    ions_per_pore = monomer.no_ions*mon_per_layer*args.layers
    for p in range(4):  # ions too
        for l in range(args.buffer, args.layers - args.buffer):
            for m in range(mon_per_layer):
                for a in range(monomer.no_ions):
                    ions.append(nmon*monomer.natoms + p*ions_per_pore + l*mon_per_layer*monomer.no_ions + m*monomer.no_ions + a)

    ref = [a.index for a in t.topology.atoms if a.name in args.ref]
    center_spline = place_solutes.trace_pores(t.xyz[0, ref, :], t.unitcell_vectors[0, :, :], args.layers)

    # Use the spline to find limits for where water is allowed in the z-direction -- only keep water inside membrane
    upper, bottom = [0, 0]
    for i in range(4):
        bottom += center_spline[i*args.layers + args.buffer, 2]
        upper += center_spline[(i + 1)*args.layers - args.buffer, 2]
    upper /= 4
    bottom /= 4

    # indices of all waters between top and bottom limits
    water = []
    for a in t.topology.atoms:
        if a.residue.name == 'HOH' and a.name == 'O':
            if bottom < t.xyz[0, a.index, 2] < upper:
                water.append(a.index)
                water.append(a.index + 1)
                water.append(a.index + 2)
    # water = [[a.index, a.index + 1, a.index + 2] for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O' and bottom < t.xyz[0, a.index, 2] < top]

    box_gromacs = [full_box[0, 0, 0], full_box[0, 1, 1], full_box[0, 2, 2], full_box[0, 0, 1], full_box[0, 2, 0],
                   full_box[0, 1, 0], full_box[0, 0, 2], full_box[0, 1, 2], full_box[0, 2, 0]]

    # Now let's work on the topology.
    # Since I'm going to place virtual sites using some, not all, of the monomers and the virtual sites may be defined
    # based on two different monomers, I will need to write the topology as a large monomer assembly
    file_rw.write_assembly(args.build_monomer, 'no', 'vsites.itp', (args.layers - 2*args.buffer)*4*mon_per_layer)
    topology = top.Top('vsites.itp')

    # I will use the benzene rings to build the virtual sites. We need to find all the benzene rings at the interface
    mon_atoms_per_pore -= (2*args.buffer*mon_per_layer*monomer.natoms)  # subtract buffered out monomers
    cndx = monomer.get_index('C')  # index of 'C' with respect to a single monomer
    c3ndx = monomer.get_index('C3')  # index of 'C3' with respect to a single monomer
    c_index = []
    c3_index = []
    for p in range(4):
        for m in range(mon_per_layer):
            c_index.append(p*mon_atoms_per_pore + m*monomer.natoms + cndx)  # bottom of pore
            # c3_index.append(p*mon_atoms_per_pore + m*monomer.natoms + c3ndx)  # bottom of pore
            # c3_index.append(p*mon_atoms_per_pore + (args.layers - args.buffer - 2)*mon_per_layer*monomer.natoms + m*monomer.natoms + c3ndx)  # top of pore

    for p in range(4):
        for m in range(mon_per_layer):
            c_index.append(p*mon_atoms_per_pore + (args.layers - 2*args.buffer - 1)*mon_per_layer*monomer.natoms + m*monomer.natoms + cndx)  # top of pore

    c3_index = [x + (c3ndx - cndx) for x in c_index]
    nreal = topology.natoms

    # add in new virtual atoms to topology [ atoms ] section
    nvsites = len(c_index)
    nbonds = nvsites // 2
    for i in range(nvsites//2):
        topology.insert_atom('v', charge=args.charge)

    for i in range(nvsites//2, nvsites):
        topology.insert_atom('v', charge=-1*args.charge)

    vndx = topology.find_indices('v')
    vndx = [i + 1 for i in vndx]

    # add some virtual bonds
    for i in range(nbonds):
        topology.insert_bond(vndx[i], vndx[i + nbonds])

    # add the actual virtual site definition to the .gro file
    vsites = np.zeros([nvsites, 3])

    stitched = t.xyz[0, mon, :]
    stitched = np.concatenate((stitched, vsites))
    stitched = np.concatenate((stitched, t.xyz[0, ions, :]))
    stitched = np.concatenate((stitched, t.xyz[0, water, :]))

    ids += ['v']*len(c_index)
    res += ['v']*len(c_index)
    res += monomer.ions*((args.layers-2*args.buffer)*mon_per_layer*4)
    ids += monomer.ions*((args.layers-2*args.buffer)*mon_per_layer*4)
    ids += ['OW', 'HW1', 'HW2']*(len(water)//3)
    res += ['SOL']*len(water)

    shrink = 3.5

    stitched[:, 2] -= (shrink/2)
    box_gromacs[2] -= shrink

    file_rw.write_gro_pos(stitched, '%s' % args.out, box=box_gromacs, ids=ids, res=res)

    for i in range(nvsites):
        nreal += 1
        topology.insert_vsite2(nreal, c_index[i], c3_index[i], 0.5)

    topology.write_top('test.itp')
