#!/usr/bin/python

"""
    The purpose of this script is to edit the topology file of a system containing molecules which have a benzene ring
    in order to create your choice of two things:
    (1) An artificial dipole which will act as electron clouds participating in a pi bond. The dipole is
        created by centering two virtual sites above and below the plane of the benzene ring and assigning them
        appropriate charges values.
    (2) Add position restraints with a given force constant to chosen atoms w.r.t. to a specified axis or axes
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

    parser.add_argument('-g', '--gro', default='initial.gro', type=str, help='Coordinate file')
    parser.add_argument('-o', '--out', default='dipole.itp', type=str, help='Name of output topology file')
    parser.add_argument('-r', '--restraints', default='on', help='Put "on" if you want position restraint on atoms')
    parser.add_argument('-D', '--dipoles', default='off', help='Put "on" if you want to create dipoles')
    parser.add_argument('-w', '--write_gro', default='off', help='Put "on" if you want to create a .gro with dipoles '
                                                                 'placed in the right spots')
    parser.add_argument('-f', '--f_const', default=1000, type=int, help='Force constant')
    parser.add_argument('-a', '--atoms', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Name of carbons in ring')
    parser.add_argument('-d', '--distance', default=0.1, help='Distance to offset dipole from ring (Angstroms)')
    parser.add_argument('-m', '--monomer', default='NAcarb11V_dummy', help='Which monomer topology is being used')
    parser.add_argument('-t', '--toplines', default=2, help='Number of lines at the top of the .gro file to ignore')
    parser.add_argument('-v', '--valence', default=1, help = 'Valence of counterion')
    parser.add_argument('-c', '--charge', default=10, help= 'Charge on dipoles')

    parser.add_argument('-A', '--axis', default='xy', help='Axis to restrain along with position restraints')

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
    ids = np.zeros([n_atoms], dtype=object)
    res = np.zeros([n_atoms], dtype=object)

    count = 0
    for i in range(top_lines, n_atoms + top_lines):
        all_coords[0, i - top_lines] = float(file[i][20:28])
        all_coords[1, i - top_lines] = float(file[i][28:36])
        all_coords[2, i - top_lines] = float(file[i][36:44])
        ids[i - top_lines] = str.strip(file[i][10:15])
        res[i - top_lines] = str.strip(file[i][5:10])
        if str.strip(file[i][10:15]) in atoms:
            coords[0, count] = float(file[i][20:28])
            coords[1, count] = float(file[i][28:36])
            coords[2, count] = float(file[i][36:44])
            count += 1


    return coords, n_atoms, all_coords, ids, res


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


def write_gro_bak(top_poles, bot_poles, b, rings):
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
        f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(count, 'PI', 'PI', count1, top_poles[0, i],
                                                    top_poles[1, i], top_poles[2, i]) + "\n")

        count += 1  # treat all virtual sites as separate residues
        count1 += 1
        f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(count, 'PI', 'NA', count1, bot_poles[0, i],
                                            bot_poles[1, i], bot_poles[2, i]) + "\n")
        count += 1
        count1 += 1
    f.write(box)
    f.close()


def write_gro(top_poles, bot_poles, b, rings, toplines, valence):
    """
    :param top_poles: The coordinates of the dipoles on the top side of the ring, numpy array [xyz coords, no rings]
    :param bot_poles: The coordinates of the dipoles on the bottom side of the ring, numpy array [xyz coords, no rings]
    :param b: A list containing each line of the original .gro file
    :param rings: the number of benzene rings
    :param toplines: number of lines at top of .gro file to ignore
    :param valence: the valence of the counter-ion used. Needed to get the counting right
    :return: Writes .gro file with new dipole locations
    """

    atoms = len(b) - 1 - toplines
    count = rings + 1  # count for the number of residues (first column of numbers in .gro file)
    count1 = atoms - (1/valence) * rings + 1  # count for total atoms (fourth column in .gro file

    for i in range(toplines):  # delete all irrelevant lines
        del b[0]

    b.insert(0, '%d\n' % (atoms + rings * 2))  # write new lines at the top of the file
    b.insert(0, 'This is a .gro file\n')

    toplines = 2

    for i in range(toplines, atoms + toplines):
        b[i] = b[i][0:5].replace(b[i][0:5], '    1') + b[i][5:len(b[i])]

    for i in range(rings):
        b.insert(count1 + toplines - 1, '{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(1, 'HII', 'PI', count1, top_poles[0, i],
                                                    top_poles[1, i], top_poles[2, i]) + "\n")

        count1 += 1
        b.insert(count1 + toplines - 1, '{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(1, 'HII', 'PI', count1, bot_poles[0, i],
                                            bot_poles[1, i], bot_poles[2, i]) + "\n")
        count1 += 1

    index = count1 + toplines - 1
    count = 2
    for i in range((1/valence)*rings):
        b[index] = '{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n'.format(count, str.strip(b[index][5:10]), str.strip(b[index][10:15]),
                                                                       count1, float(b[index][20:28]),
                                                                           float(b[index][28:36]), float(b[index][36:44]))
        index += 1
        count += 1
        count1 += 1

    f = open('dipoles.gro', 'w')  # open a new file for writing

    for line in b:  # write everything back to dipoles.gro (except the box dimension)
        f.write(line)

    f.close()


def virtual_sites(all_coords, monomers, valence, a, b, c, funct):

    n_atoms = np.shape(all_coords)[1] - (1/valence) * monomers
    atoms_per_molecule = n_atoms / monomers  # subtract valence to exclude counterion

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


def exclusions(coord_file, monomers, valence, toplines, atoms, n_atoms, vsites):
    """
    :param coord_file: the original .gro file stored in a list
    :param monomers: number of monomers
    :param valence: charge on ions
    :param atoms: the names of the atoms which should be excluded
    :param n_atoms: the number of atoms total
    :return: a list of exclusions
    """
    n_atoms = np.shape(all_coords)[1] - (1/valence) * monomers
    atoms_per_molecule = n_atoms / monomers  # subtract valence to exclude counterion (includes new PI atoms)

    n_excluded = len(atoms) + 1  # number of atoms excluded. +1 because it will be excluded from it complementary vsite

    n_exclusions = monomers * 2  # number of exclusions to be specified (one for each vsite)
    exclusions = np.zeros([n_excluded + 1, n_exclusions])  # the first entry is the virtual site itself

    for i in range(np.shape(vsites)[1]):
        exclusions[0, i] = vsites[0, i]
        # add exclusion from the complementary vsite. This is pretty specific to the format and will likely need to be
        # re-written eventually
        if i % 2 == 0:
            exclusions[1, i] = vsites[0, i + 1]
        elif i % 2 == 1:
            exclusions[1, i] = vsites[0, i - 1]

    x = 0
    for i in range(monomers):
        a = 2
        for j in range(atoms_per_molecule):
            line = i*atoms_per_molecule + j + toplines
            if str.strip(coord_file[line][10:15]) in atoms:
                exclusions[a, x] = int(coord_file[line][15:20])
                exclusions[a, x + 1] = int(coord_file[line][15:20])
                a += 1
        x += 2

    return exclusions


def position_restraints(file, atoms, axis):
    """
    Restrain the selected atoms in desired directions
    :param file: a list where each entry is a line from a coordinate file (.gro)
    :param atoms: the atoms to positions restrain
    :param axis: which direction to restrain
    :return: an array of position restraints formatted for easy writing into the topology (.itp)
    """

    # define force constants in their respective directions
    fcx = 0
    fcy = 0
    fcz = 0
    if 'x' in axis:
        fcx = args.f_const  # a large enough restraint to cause a large movement penalty
    if 'y' in axis:
        fcy = args.f_const
    if 'z' in axis:
        fcz = args.f_const

    atom_numbers = []  # find the numbers of the atoms which we are restraining
    for line in file:
        if str.strip(line[10:15]) in atoms:
            atom_numbers.append(int(line[15:20]))

    restraints = np.zeros([5, len(atom_numbers)], dtype=int)  # organize them into a list which can be translated to a topology
    for i in range(len(atom_numbers)):
        restraints[:, i] = [atom_numbers[i], 1, fcx, fcy, fcz]  # See: http://www.gromacs.org/Documentation/How-tos/Position_Restraints

    return restraints

if __name__ == "__main__":

    args = initialize()

    f = open(args.gro, 'r')
    gro = []
    for line in f:
        gro.append(line)
    f.close()

    coords, atoms, all_coords = get_coordinates(gro, 2, args.atoms)[:3]

    centers = ring_center(coords, args.atoms)

    rings = np.shape(centers)[1]

    location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Location of this script

    f = open("%s/../Structure-Files/Monomer_Tops/%s" % (location, '%s.itp' % args.monomer), "r")

    a = []
    for line in f:
        a.append(line)

    f.close()

    Assembly_itp.write_file(a, 'on', args.out, rings)

    if args.dipoles == 'on':

        perp_vectors = perpendicular_vectors(coords, args.atoms)

        top_poles, bot_poles = dipole_points(centers, perp_vectors, float(args.distance)/10) #convert distance to nanometers

        valence = int(args.valence)

        if args.write_gro == 'on':
            write_gro(top_poles, bot_poles, gro, rings, 2, valence)

        # Parameters for virtual site :
        #      /i\        i, j and k are the points from the plane chosen
        #     /   \
        #    /  C  \
        #  j/______ \k
        # Using every other carbon in the benzene rings makes the above triangle equilateral which corresponds to the
        # following parameters

        a = 1.0 / 3.0  # distance along vector rij to reach the same 'height' as C
        b = 1.0 / 3.0  # distance along vector rik to reach the same 'height' as C
        c = float(args.distance)  # The distance out of the plane to place a virtual site
        funct = 4  # this specifies the 3out type of virtual site
        vsites = virtual_sites(all_coords, np.shape(centers)[1], valence, a, b, c, funct)
        excluded_atoms = ['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O', 'O1', 'O2', 'C7', 'C21', 'C35', 'H', 'H1', 'H2',
                          'H3', 'H27', 'H28', 'H52', 'H53', 'C8', 'C22', 'C36']
        exclusions = exclusions(gro, np.shape(centers)[1], valence, 2, excluded_atoms, atoms, vsites)

        f = open(args.out, 'r')
        a = []
        atoms_count = 0
        for line in f:
            a.append(line)

        f.close()
        for i in range(len(a) - 1):
            if a[i].count('[ atoms ]') == 1:
                while a[atoms_count] != '\n':
                    atoms_count += 1
                break
            atoms_count += 1

        charge = float(args.charge)
        for i in range(rings*2):
            atom_no = atoms + i + 1 - (1/valence)*rings
            a.insert(atoms_count, '{:5d}{:>5s}{:6d}{:>6s}{:>6s}{:7d}{:5s}{:>1.6f}{:5s}{:2.5f}'.format(atom_no,
                                    'PI', 1, 'HII', 'PI', atom_no,'', (-1)**i * charge,'', 0.00000) + "\n")
            atoms_count += 1

        f = open(args.out, 'w')

        for line in a:
            f.write(line)

        f.write("\n[ virtual_sites3 ]\n")
        f.write("; Site   from                  funct         a             b             c" + "\n")
        for i in range(rings*2):
            f.write('{:<8d}{:<8d}{:<8d}{:<8d}{:<8d}{:<1.9f}{:5s}{:<1.9f}{:5s}{:<1.9f}'.format(int(vsites[0, i]), int(vsites[1, i]),
                                                                    int(vsites[2, i]), int(vsites[3, i]), int(vsites[4, i]),
                                                                    vsites[5, i],'', vsites[6, i],'', vsites[7, i]) + "\n")

        f.write("\n[ exclusions ]\n")
        for i in range(rings*2):
            for j in range(np.shape(exclusions)[0]):
                f.write('{:<8d}'.format(int(exclusions[j, i])))
            f.write("\n")

        f.close()

    if args.restraints == 'on':

        restraints = position_restraints(gro, args.atoms, '%s' % args.axis)

        f = open(args.out, 'a')  # 'a' means append

        f.write("\n[ position_restraints ]\n")
        for i in range(restraints.shape[1]):
            f.write('{:6d}{:5d}{:7d}{:7d}{:7d}'.format(restraints[0, i], restraints[1, i], restraints[2, i],
                                                       restraints[3, i], restraints[4, i]) + "\n")

        f.close()

    print 'dipole.itp file written :)'