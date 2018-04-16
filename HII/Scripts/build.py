#! /usr/bin/env python

import numpy as np
import argparse
from llclib import file_rw
from llclib import transform
import mdtraj as md
import lc_class
import os


def initialize():

    parser = argparse.ArgumentParser(description='Build HII LLC unit cell')

    parser.add_argument('-b', '--build_mon', default='NAcarb11V.gro', type=str, help='Name of class of monomer using to build with')
    parser.add_argument('-o', '--out', default='initial.gro', help='name of output file')
    parser.add_argument('-l', '--layers', default=20, type=int, help='Number of Layers')
    parser.add_argument('-m', '--monomers', default=5, type=int, help='Monomers per layer')
    parser.add_argument('-r', '--radius', default=6, type=float, help='Initial Pore Radius (Angstroms)')
    parser.add_argument('-p', '--p2p', default=45, type=float, help='Initial Pore to Pore Distance')
    parser.add_argument('-n', '--nopores', default=4, type=int, help='Number of Pores')
    parser.add_argument('-d', '--dbwl', default=3.7, type=float, help='Distance between layers')
    parser.add_argument('-s', '--layer_distribution', default='uniform', help='The distribution of monomers per layer')
    parser.add_argument('-a', '--alt_1', default=6, type=int, help='Monomers per layer for the first type of alternating layer')
    parser.add_argument('-A', '--alt_2', default=8, type=int, help='Monomers per layer for the second type of alternating layer')
    parser.add_argument('-t', '--tilt', default=0, type=float, help='Monomer tilt angle')
    parser.add_argument('-H', '--helix', help="Specify this flag if you want to build in a helical configuration",
                        action="store_true")
    parser.add_argument('-O', '--offset', help="Specify this flag to build the system in an offset configuration",
                        action="store_true")
    parser.add_argument('--rot', default=45, type=float, help="Rotate pores by this amount (degrees)")
    parser.add_argument('--offset_angle', default=0, type=float)
    parser.add_argument('-box', '--box_lengths', nargs='+', type=float, help='box vector lengths '
                                                                                                '[x y z] ')
    parser.add_argument('-angles', '--angles', nargs='+', default=[90, 90, 60], type=float, help='angles between box'
                                                                                                  'box vectors')
    parser.add_argument('-rd', '--radial_displacement', type=float, help='Shift pore center by this value for every'
                                                                         'other layer')
    parser.add_argument('-rdtheta', default=0, type=float, help='When radially displaced, the angle with respect'
                        'to the x-axis defining the direction in which to shift the pore center')
    parser.add_argument('-ad', '--angularly_displaced', action="store_true", help='rotate phenyl groups'
                        'with respect to those in adjacent layers')
    parser.add_argument('-columns', action="store_true", help='build column-wise')
    parser.add_argument('-L', '--correlation_length', default=10, type=float, help='Desired correlation length')
    parser.add_argument('-Lvar', default=0.1, type=float, help='variance in z position of monomer heads')

    args = parser.parse_args()

    return args


def z_correlation(z, L, v=0.1):
    """
    Calculate where to place monomers on the z-axis so that a given correlation length is obtained
    :param z: mean z-positions where monomers will be placed with gaussian probability np.array([n_layers])
    :param L: desired correlation length [float]
    :param v: variance in z position of monomer head groups
    :return: locations [np.array[nlayers])
    """

    n = z.shape[0]
    cov = np.zeros([n, n])  # initialize covariance matrix

    decay = v*np.exp(-z / L)  # decay of covariance
    # decay[1:] += np.exp(-z[::-1][:-1]/L) # for periodicity (?)

    for i in range(z.shape[0]):
        cov[i, i:] = decay[:(n - i)]
        cov[i:, i] = decay[:(n - i)]

    locations = np.random.multivariate_normal(z, cov)

    return locations


if __name__ == "__main__":

    args = initialize()

    props = lc_class.LC('%s' % args.build_mon)

    no_monomers = args.monomers  # number of monomers packed per layer around a pore
    pore_radius = args.radius  # Radius of pore (unsure of units right now)
    no_pores = args.nopores  # number of pores to be simulated
    p2p = args.p2p / 10  # distance between pores (units tbd)
    no_layers = args.layers  # Number of layers in a pore
    dist = args.dbwl  # distance between layers (units tbd)
    nmon = no_pores*no_layers*no_monomers

    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    t = md.load("%s/../top/HII_Monomer_Configurations/%s" % (location, args.build_mon))
    pos = t.xyz[0, :, :]

    natoms = nmon * pos.shape[0]
    res = [a.residue.name for a in t.topology.atoms]
    ids = [a.name for a in t.topology.atoms]

    if args.angularly_displaced:
        phenyl_carbons = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
        phenyl_indices = [a.index for a in t.topology.atoms if a.name in phenyl_carbons]

    if args.box_lengths:
        a, b, c = args.box_lengths
    else:
        a, b, c = [2*p2p, 2*p2p, args.dbwl*args.layers/10]  # logical choices for symmetry

    alpha, beta, gamma = [angle * (np.pi / 180) for angle in args.angles]

    V = a * b * c * np.sqrt(
        1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.sin(
            gamma)) # volume of unitcell

    A = np.array([a, 0, 0])  # vector in x direction
    B = np.array([b*np.cos(gamma), b*np.sin(gamma), 0])  # vector in y direction
    C = np.array([c*np.cos(beta), c*((np.cos(alpha) - np.cos(gamma)*np.cos(beta))/np.sin(gamma)), V / (a*b*np.sin(gamma))])  # vector in z direction

    box = [A[0], B[1], C[2], A[1], A[2], B[0], B[2], C[0], C[1]]  # gromacs box format

    # rotate monomer so plane of aromatic head groups is coplanar with xy plane

    plane_atoms = np.zeros([3, 3])
    for i in range(plane_atoms.shape[0]):
        plane_atoms[i, :] = pos[props.plane_indices[i], :]

    R = transform.rotateplane(plane_atoms, angle=args.tilt)  # generate rotation matrix

    b = np.ones([1])
    for i in range(pos.shape[0]):
        coord = np.concatenate((pos[i, :], b))
        x = np.dot(R, coord)
        pos[i, :] = x[:3]

    # translate molecule to origin
    pos = transform.translate(pos, pos[props.ref_atom_index, :], [0, 0, 0])

    # align monomer with x-axis
    v = np.array([pos[props.lineatoms[0], :2] - pos[props.lineatoms[1], :2]])
    angle = np.arctan(v[0, 1]/v[0, 0])
    pos = transform.rotate_coords_z(pos, -angle*180/np.pi)
    v = np.array([pos[props.lineatoms[0], :2] - pos[props.lineatoms[1], :2]])

    # calculate locations of pore centers for monoclinic setup
    pore_centers = np.zeros([4, 2])

    # theta = np.pi * args.angles[2] / 180
    cell_theta = np.pi / 3

    pore_centers[0, :] = [p2p/2, p2p/2]
    pore_centers[1, :] = [p2p/2, 3*p2p/2]
    pore_centers[2, :] = [3*p2p/2, 3*p2p/2]
    pore_centers[3, :] = [3*p2p/2, p2p/2]

    pore_centers[..., 1] = pore_centers[..., 1] * np.sin(cell_theta)
    pore_centers[..., 0] = pore_centers[..., 0] + pore_centers[..., 1]*np.cos(cell_theta)

    # place monomers
    assembly = np.zeros([pos.shape[0] * args.layers * args.monomers * args.nopores, 3])

    wedge_angle = 2 * np.pi / args.monomers
    r = args.radius / 10

    if args.columns:

        mean_z = np.linspace(0, (args.layers - 1) * args.dbwl, args.layers) / 10  # mean z-positions of monomer heads
        mean_z[::2] -= 0.15
        rotated_carboxylate = np.copy(pos)

        rotated_carboxylate[props.carboxylate_indices, :] = transform.rotate_coords_x(
            rotated_carboxylate[props.carboxylate_indices, :], 90)

        for p in range(no_pores):
            px = pore_centers[p, 0]
            py = pore_centers[p, 1]
            for m in range(no_monomers):
                theta = m * wedge_angle + args.rot * np.pi / 180
                #offset_theta = theta + wedge_angle/4
                antiparallel = transform.rotate_coords_z(pos, theta*180/np.pi)
                parallel = transform.rotate_coords_z(rotated_carboxylate, theta*180/np.pi)
                shift = np.random.uniform(-args.dbwl/20, args.dbwl/20)
                #z = z_correlation(mean_z, args.correlation_length, args.Lvar)
                # z = mean_z + shift
                z = mean_z
                for l in range(z.shape[0]):
                    ndx = p*args.layers*args.monomers*pos.shape[0] + m*args.layers*pos.shape[0] + l*pos.shape[0]
                    if l % 2 == 0:
                        assembly[ndx:(ndx + pos.shape[0]), :] = transform.translate(parallel,
                                    pos[props.ref_atom_index, :], [px + r*np.cos(theta), py + r*np.sin(theta), z[l]])
                    else:
                        assembly[ndx:(ndx + pos.shape[0]), :] = transform.translate(antiparallel,
                                    pos[props.ref_atom_index, :], [px + r*np.cos(theta), py + r*np.sin(theta), z[l]])

    else:

        # calculate locations of layers in z-direction
        z = np.linspace(0, args.dbwl*(args.layers - 1)/10, args.layers)

        for p in range(args.nopores):
            for l in range(args.layers):
                for m in range(args.monomers):
                    if args.offset and l % 2 == 0:
                        theta = m*wedge_angle + wedge_angle/2 + args.rot*np.pi/180
                    else:
                        theta = m*wedge_angle + args.rot*np.pi/180
                    if args.radial_displacement and l % 2 == 0:
                        px = pore_centers[p, 0] + args.radial_displacement*np.cos(np.pi*args.rdtheta/180)
                        py = pore_centers[p, 1] + args.radial_displacement*np.sin(np.pi*args.rdtheta/180)
                    else:
                        px = pore_centers[p, 0]
                        py = pore_centers[p, 1]
                    ndx = p*args.monomers*args.layers*pos.shape[0] + l*args.monomers*pos.shape[0] + m*pos.shape[0]
                    # rotate monomer about origin by ref-atom
                    rotated = transform.rotate_coords_z(pos, theta*180/np.pi)
                    if args.angularly_displaced and l % 2 == 0:
                        phenyl_centroid = np.mean(rotated[phenyl_indices, :], axis=0)  # center of phenyl ring
                        trans = transform.translate(rotated, phenyl_centroid, [0, 0, 0])  # translate center to origin
                        rot = transform.rotate_coords_z(trans, 30)  # rotate monomer about z-axis by 30 degrees
                        rotated = transform.translate(rot, [0, 0, 0], phenyl_centroid)  # translate centroid back
                    # translate monomer away from center by the pore radius
                    rotated[:, :2] += [r*np.cos(theta), r*np.sin(theta)]
                    # translate monomer to appropriate pore center and z height
                    assembly[ndx:(ndx + pos.shape[0]), :] = transform.translate(rotated, pos[props.ref_atom_index, :],
                            [px, py, z[l]])

    # reorder so that monomer and ion residues are separate
    ions = []
    for n in range(props.no_ions):
        for i in range(args.monomers*args.layers*args.nopores):
            ions.append(props.ion_indices[n]*(i + 1) + props.no_ions*i)

    monomer = [i for i in range(natoms) if i not in ions]

    ordered = monomer + ions

    assembly = assembly[ordered, :]

    res = np.array(res*nmon)
    ids = np.array(ids*nmon)

    res = res[ordered]

    ids = ids[ordered]

    file_rw.write_gro_pos(assembly, '%s' % args.out, res=list(res), ids=list(ids), box=box)
