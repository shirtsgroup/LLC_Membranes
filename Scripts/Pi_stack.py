#!/usr/bin/bash

"""
    The purpose of this script is to edit the topology file of a system containing molecules which have a benzene ring
    in order to create an artificial dipole which will as electron clouds participating in a pi bond. The dipole is
    created by centering two virtual sites above and below the plane of the benzene ring and assigning them appropriate
    charges values.
"""

import argparse
import Get_Positions
import numpy as np
import math
import warnings
import Assembly_itp
import subprocess
import os
import Assembly_itp


def initialize():

    parser = argparse.ArgumentParser(description='Duplicate points periodically in the x-y directions')

    parser.add_argument('-f', '--file', default='NaPore.top', help='File to replicate periodically')
    parser.add_argument('-c', '--coord_file', default='test.gro', help='Coordinate file')
    parser.add_argument('-o', '--output', default='NaPore_Pi.top', help='Name of output file')
    parser.add_argument('-a', '--atoms', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Name of carbons in ring')
    parser.add_argument('-d', '--distance', default=0.1, help='Distance to offset dipole from ring (Angstroms)')
    parser.add_argument('-m', '--monomer', default='monomer2', help='Which monomer topology is being used')
    parser.add_argument('-t', '--toplines', default=2, help='Number of lines at the top of the .gro file to ignore')
    parser.add_argument('-v', '--valence', default=1, help = 'Valence of counterion')

    args = parser.parse_args()

    return args

warnings.filterwarnings("error")  # This makes it so numpy warnings are treated as real errors


def get_coordinates(file, top_lines, atoms):
    """
    Get the coordinates of all atoms and atoms of interest from a .gro file
    :param file: a list with each entry containing a line from the .gro file
    :param top_lines: the number of lines at the top of the .gro file to skip
    :param atoms: A list containing residue names of atoms in benzene rings, user input
    :return: Coordinates of all benzene carbons (numpy array), the total number of atoms in the system (int),
    the coordinates of all atoms (numpy array)
    """

    n_atoms = len(file) - top_lines - 1  # subtract 1 more for the box at the bottom

    count = 0
    for i in range(top_lines, n_atoms + top_lines):
        if str.strip(file[i][10:15]) in atoms:
            count += 1

    coords = np.zeros([3, count])
    all_coords = np.zeros([3, n_atoms])

    count = 0
    for i in range(top_lines, n_atoms + top_lines):
        all_coords[0, i - top_lines] = float(file[i][20:28])
        all_coords[1, i - top_lines] = float(file[i][28:36])
        all_coords[2, i - top_lines] = float(file[i][36:44])
        if str.strip(file[i][10:15]) in atoms:
            coords[0, count] = float(file[i][20:28])
            coords[1, count] = float(file[i][28:36])
            coords[2, count] = float(file[i][36:44])
            count += 1

    return coords, n_atoms, all_coords


def ring_center(coords, atoms):
    """
    :param coords: coordinates of all benzene atoms, numpy array [xyz positions, number of benzene atoms]
    :param atoms: A list containing residue names of atoms in benzene rings, user input
    :return: The center coordinates of each benzene ring, numpy array
        NOTE: This is actually unnecessary to construct a virtual site in gromacs. Setting up the virtual site properly will
    take care of the placement of the virtual site. You can initially place the virtual site 'atoms' anywhere and their
    positions will be corrected within the first few steps of the simulation
    """
    no_atoms = len(atoms)
    rings = np.shape(coords)[1] / no_atoms
    centers = np.zeros([3, rings])

    for i in range(rings):
        sum = np.array([0.0, 0.0, 0.0])
        for j in range(no_atoms):
            sum += coords[:, i*no_atoms + j]
        centers[:, i] = sum / float(no_atoms)

    return centers


def perpendicular_vectors(coords, atoms):
    """
    :param coords: The coordinates of all benzene carbons in the system, numpy array [xyz positions, no benzene atoms]
    :param atoms: A list containing residue names of atoms in benzene rings, user input
    :return: Perpendicular vectors to the plane of atoms in the benzene rings
        NOTE: This is actually unnecessary to construct a virtual site in gromacs. Setting up the virtual site properly will
    take care of the placement of the virtual site. You can initially place the virtual site 'atoms' anywhere and their
    positions will be corrected within the first few steps of the simulation
    """
    no_atoms = len(atoms)
    rings = np.shape(coords)[1] / no_atoms
    perp_vectors = np.zeros([3, rings])

    for i in range(rings):
        pts = np.zeros([3, 3])  # need 3 points from the benzene ring plane
        for j in range(no_atoms):
            if j % 2 == 0:  # only need 3 atoms so let's take every other - this is important in the vsite construction
                pts[:, j / 2] = coords[:, i*no_atoms + j]
        v1 = pts[:, 0] - pts[:, 1]
        v2 = pts[:, 0] - pts[:, 2]
        perp_vectors[:, i] = np.cross(v1, v2)

    return perp_vectors


def dipole_points(centers, perp_vectors, distance):
    """
    This function sets the location of the virtual sites with respect to the plane formed by chosen atoms
    :param centers: coordinates of the centers of the benzene rings, numpy array [xyz position, number of rings]
    :param perp_vectors: vectors perpendicular to each plans, numpy array [xyz vector, number of rings]
    :param distance: the distance to offset virtual sites from the plane of atoms, user input
    :return: The coordinates of the locations of all dipoles, numpy array
        NOTE: This is actually unnecessary to construct a virtual site in gromacs. Setting up the virtual site properly will
    take care of the placement of the virtual site. You can initially place the virtual site 'atoms' anywhere and their
    positions will be corrected within the first few steps of the simulation
    """
    rings = np.shape(centers)[1]
    top_pole = np.zeros([3, rings])
    bot_pole = np.zeros([3, rings])

    for i in range(rings):
        perp_vect = perp_vectors[:, i]
        center = centers[:, i]
        try:  # this is the reason warnings are filtered. Numpy encounters a divide by zero if the plane is perfectly
              # in the x-y plane as it is when making an initial configuration
            norm = perp_vect / math.sqrt(np.dot(perp_vect, perp_vect))
        except RuntimeWarning:
            norm = np.array([0.0, 0.0, 1.0])
        if norm[2] < 0:
            top_pole[:, i] = center - distance*norm
            bot_pole[:, i] = center + distance*norm
        else:
            top_pole[:, i] = center + distance*norm
            bot_pole[:, i] = center - distance*norm

    return top_pole, bot_pole


def write_gro(top_poles, bot_poles, b, rings):
    """
    :param top_poles: The coordinates of the dipoles on the top side of the ring, numpy array [xyz coords, no rings]
    :param bot_poles: The coordinates of the dipoles on the bottom side of the ring, numpy array [xyz coords, no rings]
    :param b: A list containing each line of the original .gro file
    :param rings: the number of benzene rings
    :return: Writes .gro file with new dipole locations
    """

    f = open('dipoles.gro', 'w')  # open a new file for writing

    box = b[len(b) - 1]  # extract the box dimensions
    del b[len(b) - 1]  # delete the box dimensions from the file

    for line in b:  # write everything back to dipoles.gro (except the box dimension)
        f.write(line)

    count = rings * 2 + 1  # count for the number of residues (first column of numbers in .gro file)
    count1 = atoms + 1  # count for total atoms (fourth column in .gro file
    for i in range(rings):
        f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(count, 'NA', 'NA', count1, top_poles[0, i],
                                                    top_poles[1, i], top_poles[2, i]) + "\n")

        count += 1  # treat all virtual sites as separate residues
        count1 += 1
        f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(count, 'NA', 'NA', count1, bot_poles[0, i],
                                            bot_poles[1, i], bot_poles[2, i]) + "\n")
        count += 1
        count1 += 1
    f.write(box)
    f.close()


def virtual_sites(all_coords, monomers, valence, a, b, c, funct):

    n_atoms = np.shape(all_coords)[1]
    atoms_per_molecule = n_atoms / monomers - valence  # subtract valence to exclude counterion
    n_vsites = monomers * 2  # number of virtual sites need to create a dipole (2 per ring)
    vsites = np.zeros([8, n_vsites])
    vsites[4, :] = funct
    vsites[5, :] = a  # all a and b values are the same
    vsites[6, :] = b

    for i in range(monomers):
        for j in range(2):
            vsites[0, i * 2 + j] = n_atoms + i * 2 + j + 1  # new atoms were placed at the end of the .gro file
            vsites[1:4, i * 2 + j] = [i * atoms_per_molecule + 1, i * atoms_per_molecule + 3,
                              i * atoms_per_molecule + 5]
            vsites[7, i * 2 + j] = (-1)**j * c

    return vsites


if __name__ == "__main__":

    args = initialize()

    f = open('%s' % args.file, 'r')
    a = []
    for line in f:
        a.append(line)
    f.close()

    f = open('%s' % args.coord_file, 'r')
    b = []
    for line in f:
        b.append(line)
    f.close()

    coords, atoms, all_coords = get_coordinates(b, 2, args.atoms)

    centers = ring_center(coords, args.atoms)

    perp_vectors = perpendicular_vectors(coords, args.atoms)

    top_poles, bot_poles = dipole_points(centers, perp_vectors, float(args.distance))

    rings = np.shape(centers)[1]

    write_gro(top_poles, bot_poles, b, rings)

    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Location of this script

    f = open("%s/../Structure-Files/Monomer_Tops/%s" % (location, '%s.itp' % args.monomer), "r")

    a = []
    for line in f:
        a.append(line)

    f.close()

    Assembly_itp.write_file(a, 'off', 'dipole.itp', rings)
    print 'dipole.itp file written'

    # Parameters for virtual site :
    #      /i\        i, j and k are the points from the plane chosen
    #     /   \
    #    /  C  \
    #  j/______ \k
    # Using every other carbon in the benzene rings makes the above triangle equilateral which corresponds to the
    # following parameters
    a = 2 / math.sqrt(3)  # distance along vector rij to reach the same 'height' as C
    b = 2 / math.sqrt(3)  # distance along vector rik to reach the same 'height' as C
    c = float(args.distance)  # The distance out of the plane to place a virtual site
    funct = 4  # this specifies the 3out type of virtual site

    vsites = virtual_sites(all_coords, np.shape(centers)[1], int(args.valence), a, b, c, funct)

    f = open('dipole.itp', 'a')

    f.write("[ virtual_sites3 ]" + "\n")
    f.write("; Site   from                  funct         a             b             c" + "\n")
    for i in range(rings*2):
        f.write('{:<8d}{:<8d}{:<8d}{:<8d}{:<8d}{:<1.9f}{:5s}{:<1.9f}{:5s}{:<1.9f}'.format(int(vsites[0, i]), int(vsites[1, i]),
                                                                int(vsites[2, i]), int(vsites[3, i]), int(vsites[4, i]),
                                                                vsites[5, i],'', vsites[6, i],'', vsites[7, i]) + "\n")

    f.close()