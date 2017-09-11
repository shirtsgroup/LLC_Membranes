#! /usr/bin/env python

"""
Rotate monomer tails so they form a specified angle with the xy plane
"""

import argparse
import numpy as np
import mdtraj as md
import reposition
import tilt


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the tilt angle of alkyl tails')

    parser.add_argument('-m', '--mon', default='tilted.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-o', '--output', default='tilt.gro', type=str, help ='Name of output file')
    parser.add_argument('-b', '--bond_axes', nargs='+', default=["C4-O", "C3-O1", "C2-O2"], help='Bonds around which to rotate')
    parser.add_argument('-a', '--angles', nargs='+', default=[40, 40, 40], type=float, help='Trajectory file (.xtc or .trr)')
    parser.add_argument('--gen_tails', help='Convenience flag for generating an index group containing the tail'
                                            'indices. This only works with the monomer which this program was written'
                                            'for', action="store_true")
    parser.add_argument('-x', '--index', type = str, default='index.ndx', help='Index file containing groups for each '
                        'tail. Each group should have a header line of the format [ groupname ]. The following line '
                        'should list the atoms in the order of their connectivity Each group should be separated by a '
                                                                               'blank line')

    args = parser.parse_args()

    return args


def atoms(bond_list, top):
    """
    :param bonds: A list of names of atoms in bonded pairs
    :return: A list with just atom names
    """

    nbonds = len(bond_list)  # number of bonds

    atoms = []  # reformat input from ["C-C"] to ["C", "C"]
    for i in range(nbonds):
        a = bond_list[i].split('-')
        for j in range(len(a)):
            atoms.append(a[j])

    indices = np.zeros([nbonds*2], dtype=int)
    for atom in top:
        if atom.name in atoms:
            ndx = atoms.index(atom.name)
            indices[ndx] = atom.index

    indices = np.reshape(indices, (nbonds, 2))

    return indices


def tail_rotate(axis, tail, angle):

    # technique taken from here: http://paulbourke.net/geometry/rotate/
    angle *= (np.pi / 180)
    natoms = tail.xyz.shape[1]

    # Step 1: Translate point on vector to origin

    T = reposition.translate(axis[0, :])
    Tinv = np.copy(T)
    Tinv[:3, 3] = -T[:3, 3]

    # Step 2: Rotate space about the x axis so that the rotation axis lies in the xz plane

    u = axis[0, :] - axis[1, :]
    U = u / np.linalg.norm(u)
    a = U[0]
    b = U[1]
    c = U[2]
    d = np.sqrt(b**2 + c**2)
    Rx = np.zeros([4, 4])
    Rx[0, 0] = 1
    Rx[3, 3] = 1
    Rx[1, 1] = c/d
    Rx[1, 2] = -b/d
    Rx[2, 1] = b/d
    Rx[2, 2] = c/d
    Rxinv = np.copy(Rx)
    Rxinv[1, 2] = b/d
    Rxinv[2, 1] = -b/d

    # Step 3: Rotate space about the y axis so that the rotation axis lies along the positive z axis

    Ry = np.zeros([4, 4])
    Ry[0, 0] = d
    Ry[1, 1] = 1
    Ry[2, 2] = d
    Ry[3, 3] = 1
    Ry[0, 2] = -a
    Ry[2, 0] = a
    Ryinv = np.copy(Ry)
    Ryinv[0, 2] = a
    Ryinv[2, 0] = -a

    # Step 4: Rotation about the z-axis by theta

    Rz = np.zeros([4, 4])
    Rz[0, 0] = np.cos(angle)
    Rz[0, 1] = -np.sin(angle)
    Rz[1, 0] = np.sin(angle)
    Rz[1, 1] = np.cos(angle)
    Rz[2, 2] = 1
    Rz[3, 3] = 1

    # rotate = T.dot(Rx).dot(Ry).dot(Rz).dot(Ryinv).dot(Rxinv).dot(Tinv)
    rotate = Tinv.dot(Rxinv).dot(Ryinv).dot(Rz).dot(Ry).dot(Rx).dot(T)
    # Apply all these transformations to each point
    # rotated = np.zeros([natoms, 3])
    for i in range(natoms):
        cat = np.append(tail.xyz[0, i, :], [1])
        tail.xyz[0, i, :] = np.dot(rotate, cat)[:3]

    return tail


def write_gro(full, out):

    pos = full.xyz

    ids = [a.name for a in full.topology.atoms]
    res = []
    for i in range(full.xyz.shape[1]):
        res.append(full.topology.atom(i).residue.name)

    box = full.unitcell_vectors

    atoms = pos.shape[1]
    count = 1
    with open(out, 'w') as f:
        f.write('Solvated system with water layer removed\n')
        f.write('%s\n' % atoms)
        for i in range(atoms):
            f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(count, res[i], ids[i],
                        count, pos[0, i, 0], pos[0, i, 1], pos[0, i, 2]))
            count += 1
        f.write('{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}{:10f}\n'.format(box[0, 0, 0], box[0, 1, 1], box[0, 2, 2]
                                                                    ,box[0, 0, 1], box[0, 2, 0], box[0, 1, 0]
                                                                    ,box[0, 0, 2], box[0, 1, 2], box[0, 2, 0]))


if __name__ == "__main__":

    args = initialize()

    full = md.load('%s' % args.mon)
    pos = full.xyz

    indices = atoms(args.bond_axes, full.topology.atoms)  # get indices of atoms in bond which we are rotating about

    tails = tilt.read_index(args.index)

    ntails = len(tails)
    rotated = np.zeros([ntails, len(tails[0]), 3])

    for i in range(ntails):

        keep = [a.index for a in full.topology.atoms if a.name in tails[i]]
        t = full.atom_slice(keep)

        ax = np.zeros([2, 3])  # points making up rotation axis
        ax[0, :] = full.xyz[0, indices[i][0], :]
        ax[1, :] = full.xyz[0, indices[i][1], :]

        t = tail_rotate(ax, t, args.angles[i])

        # update the positions in full.xyz

        for i in range(t.xyz.shape[1]):
            for j in range(full.xyz.shape[1]):
                if t.topology.atom(i).name == full.topology.atom(j).name:
                    full.xyz[0, j, :] = t.xyz[0, i, :]

    write_gro(full, args.output)







