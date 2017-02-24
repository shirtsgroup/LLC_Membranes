#!/usr/bin/env python

import numpy as np
import math
import os
import argparse
import LC_class

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def initialize():

    parser = argparse.ArgumentParser(description='Build LLC Structure')

    parser.add_argument('-b', '--build_mon', default='NAcarb11Vd', type=str, help='Name of class of monomer using to build with')
    parser.add_argument('-o', '--out', default='initial.gro', help='name of output file')
    parser.add_argument('-l', '--layers', default=20, type=int, help = 'Number of Layers')
    parser.add_argument('-m', '--monomers', default=6, type=int, help = 'Monomers per layer')
    parser.add_argument('-r', '--radius', default=6, type=float, help = 'Initial Pore Radius (Angstroms)')
    parser.add_argument('-p', '--p2p', default=45, type=float, help = 'Initial Pore to Pore Distance')
    parser.add_argument('-n', '--nopores', default=4, type=int, help = 'Number of Pores')
    parser.add_argument('-d', '--dbwl', default=5, type=float, help = 'Distance between layers')
    parser.add_argument('-s', '--layer_distribution', default='uniform', help = 'The distribution of monomers per layer')
    parser.add_argument('-a', '--alt_1', default=6, type=int, help='Monomers per layer for the first type of alternating layer')
    parser.add_argument('-A', '--alt_2', default=8, type=int, help='Monomers per layer for the second type of alternating layer')
    parser.add_argument('-t', '--tilt', default=0, type=float, help='Monomer tilt angle')
    parser.add_argument('-H', '--helix', help="Specify this flag if you want to build in a helical configuration",
                        action="store_true")
    parser.add_argument('-O', '--offset', help="Specify this flag to build the system in an offset configuration",
                        action="store_true")

    args = parser.parse_args()

    return args


def read_pdb_coords(file):

    a = []
    for line in file:
        a.append(line)
    f.close()

    no_atoms = 0  # number of atoms in one monomer including sodium ion
    for i in range(0, len(a)):
        no_atoms += a[i].count('ATOM')

    lines_of_text = 0  # lines of text at top of .pdb input file
    for i in range(0, len(a)):
        if a[i].count('ATOM') == 0:
            lines_of_text += 1
        if a[i].count('ATOM') == 1:
            break

    xyz = np.zeros([3, no_atoms])
    identity = np.zeros([no_atoms], dtype=object)
    for i in range(lines_of_text, lines_of_text + no_atoms):  # searches relevant lines of text in file, f, being read
        xyz[:, i - lines_of_text] = [float(a[i][26:38]), float(a[i][38:46]), float(a[i][46:54])]  # Use this to read specific entries in a text file
        identity[i - lines_of_text] = str.strip(a[i][12:16])

    return xyz, identity, no_atoms, lines_of_text


def read_gro_coords(file):

    a = []
    for line in file:
        a.append(line)
    f.close()

    lines_of_text = 2  # Hard Coded -> BAD
    no_atoms = len(a) - lines_of_text - 1  # subtract one for the bottom box vector line

    xyz = np.zeros([3, no_atoms])
    identity = np.zeros([no_atoms], dtype=object)
    for i in range(lines_of_text, lines_of_text + no_atoms):  # searches relevant lines of text in file, f, being read
        xyz[:, i - lines_of_text] = [float(a[i][20:28])*10, float(a[i][28:36])*10, float(a[i][36:44])*10]
        identity[i - lines_of_text] = str.strip(a[i][11:16])

    return xyz, identity, no_atoms, lines_of_text


def layer_dist(layers, nopores, distribution):

    layer_distribution = np.zeros([layers*nopores], dtype=int)

    if distribution == 'uniform':
        for i in range(0, len(layer_distribution)):
            layer_distribution[i] = args.monomers
    if distribution == 'alternating':
        for i in range(layer_distribution.shape[-1]):
            if i % 2 == 0:
                layer_distribution[i] = args.alt_1
            if i % 2 == 1:
                layer_distribution[i] = args.alt_2
    return layer_distribution


def slope(pt1, pt2):
    m = (pt1[1] - pt2[1])/(pt1[0] - pt2[0])  # slope
    return m


def rotateplane(plane_atoms, angle=0):
    """
    Calculate a rotation matrix to rotate a plane in 3 dimensions
    :param plane_atoms: names of atoms making up plane which is being aligned in the xy plane
    :param angle: desired angle between xy plane (optional, default = 0 i.e. in plane)
    :return:
    """
    # Arrays to hold x,y,z values of each point of interest
    plane = np.zeros([3, 3])

    count = 0
    for i in range(no_atoms):
        if identity[i] in plane_atoms:
            plane[count, :] = xyz[:, i]
            count += 1

    # vector pointing from point 1 to point 2
    v12 = plane[1, :] - plane[0, :]
    v13 = plane[2, :] - plane[0, :]

    # The cross product of v12 and v13 give a vector that is perpendicular to the plane:
    N = np.cross(v12, v13)

    N_desired = [0, math.sin(angle), math.cos(angle)]  # vector in the direction normal to our desired plane orientation

    RotationAxis = np.cross(N, N_desired)
    theta = math.acos(np.dot(N, N_desired)/(np.linalg.norm(N)*np.linalg.norm(N_desired)))  #  Rotation Angle (radians)

    L = [RotationAxis[0]/np.linalg.norm(RotationAxis), RotationAxis[1]/np.linalg.norm(RotationAxis),
                       RotationAxis[2]/np.linalg.norm(RotationAxis)]  # normalized Rotation Axis
    # ^ see: http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/

    u, v, w = L[0], L[1], L[2]

    R = np.zeros((4, 4))
    R[3, 3] = 1
    R[0, 0] = u**2 + (v**2 + w**2)*math.cos(theta)  # math.cos takes theta in radians by default
    R[0, 1] = u*v*(1 - math.cos(theta)) - w*math.sin(theta)
    R[0, 2] = u*w*(1 - math.cos(theta)) + v*math.sin(theta)
    R[1, 0] = u*v*(1 - math.cos(theta)) + w*math.sin(theta)
    R[1, 1] = v**2 + (u**2 + w**2)*math.cos(theta)
    R[1, 2] = v*w*(1 - math.cos(theta)) - u*math.sin(theta)
    R[2, 0] = u*w*(1 - math.cos(theta)) - v*math.sin(theta)  # math.cos takes theta in radians by default
    R[2, 1] = w*v*(1 - math.cos(theta)) + u*math.sin(theta)
    R[2, 2] = w**2 + (u**2 + v**2)*math.cos(theta)

    return R


def quadrant(pt, origin=[0, 0]):
    """ Find out which quadrant of the xy plane a point is sitting in
    II    |    I
          |
    -------------
          |
    III   |    IV
    :param: pt: point to be tested
    :param: origin: the location of the origin. Default is [0, 0] but can be set arbitrarily (such as a pore center)
    """
    if pt[0] > origin[0] and pt[1] > origin[1]:
        return 1
    elif pt[0] < origin[0] and pt[1] < origin[1]:
        return 3
    elif pt[0] > origin[0] and pt[1] < origin[1]:
        return 4
    elif pt[0] < origin[0] and pt[1] > origin[1]:
        return 2
    else:
        return 0  # the case where the point lies on the x or y axis


def transdir(pt, origin=[0, 0]):
    """
    figure out in which direction the coordinates will be shifted. They are always shifted away from the origin
    :param pt:
    :return:
    """
    if quadrant(pt, origin) == 1:  # e.g. in quadrant 1, the x's are shifted in the positive x and positive y directions
        vx = 1
        vy = 1
    elif quadrant(pt, origin) == 2:  # in quadrant 2, the x's are shifted negative and the y's are shifted positive
        vx = -1
        vy = 1
    elif quadrant(pt, origin) == 3:  # in quadrant 3, the x's and y's are shifted down
        vx = -1
        vy = -1
    elif quadrant(pt, origin) == 4:  # in quadrant 4, the x's are shifted positive and the y's are shifted negative
        vx = 1
        vy = -1
    # These next three conditionals are very unlikely but are included for completeness and to avoid future errors
    elif quadrant(pt, origin) == 0 and pt[0] == origin[0]:  # i.e., it lies on the y - axis
        if pt[1] > 0:  # the point is on the positive y-axis
            vx = 0  # no x-shift
            vy = 1  # shift in the positive y direction
        if pt[1] < 0:  # the point is on the negative y-axis
            vx = 0  # no x-shift
            vy = -1  # shift in the negative y direction
    elif quadrant(pt, origin) == 0 and pt[1] == origin[1]:  # i.e., it lies on the x - axis
        if pt[0] > 0:  # the point is on the positive x-axis
            vx = 1  # shift in the positive x direction
            vy = 0  # no y-shift
        if pt[0] < 0:  # the point is on the negative y-axis
            vx = -1  # shift in the negative x direction
            vy = 0  # no y-shift

    return vx, vy


def rotate(theta):
    """
    :param: angle by which to rotate the monomer
    :return: Creates a rotation matrix to rotate input vector about y-axis making a new coordinate at evenly spaced points.
    """
    Rx = np.zeros([3, 3])  # makes a 3 x 3 zero matrix
    Rx[0, 0] = math.cos(theta)
    Rx[1, 0] = math.sin(theta)
    Rx[0, 1] = -math.sin(theta)
    Rx[1, 1] = math.cos(theta)
    Rx[2, 2] = 1

    return Rx


def write_gro_bak(positions, identity, no_layers, layer_distribution, dist, no_pores, p2p, no_ions):

    f = open('%s' % args.out, 'w')

    f.write('This is a .gro file\n')
    f.write('%s\n' % sys_atoms)

    atom_count = 1
    monomer_count = 0
    for l in range(0, no_pores):  # loop to create multiple pores
        theta = 30  # angle which will be used to do hexagonal packing
        if l == 0:  # unmodified coordinates
            b = 0
            c = 0
        elif l == 1:  # move a pore directly down
            b = - 1
            c = 0
        elif l == 2:  # moves pore up and to the right
            b = math.cos(math.radians(90 - theta))
            c = -math.sin(math.radians(90 - theta))
        elif l == 3:  # moves a pore down and to the right
            b = -math.sin(math.radians(theta))
            c = -math.cos(math.radians(theta))
        for k in range(0, no_layers):
            layer_mons = layer_distribution[l*no_layers + k]
            positions_index = sorted(set(layer_distribution)).index(layer_mons)  # The index in positions which should be read from
            for j in range(0, layer_mons):  # iterates over each monomer to create coordinates
                monomer_count += 1
                for i in range(0, no_atoms - no_ions):
                    x_values.append(b*p2p + positions[positions_index][j][i][0])
                    y_values.append(c*p2p + positions[positions_index][j][i][1])
                    z_values.append(k*dist + + (dist/float(layer_mons))*j + positions[positions_index][j][i][2])
                    hundreds = int(math.floor(atom_count/100000))
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(monomer_count, name, identity[i],
                        atom_count - hundreds*100000, x_values[atom_count - 1]/10.0, y_values[atom_count - 1]/10.0,
                        z_values[atom_count - 1]/10.0) + "\n")
                    atom_count += 1

    # Ions:

    for l in range(0, no_pores):  # loop to create multiple pores
        theta = 30  # angle which will be used to do hexagonal packing
        if l == 0:  # unmodified coordinates
            b = 0
            c = 0
        elif l == 1:  # move a pore directly down
            b = - 1
            c = 0
        elif l == 2:  # moves pore up and to the right
            b = math.cos(math.radians(90 - theta))
            c = -math.sin(math.radians(90 - theta))
        elif l == 3:  # moves a pore down and to the right
            b = -math.sin(math.radians(theta))
            c = -math.cos(math.radians(theta))
        for k in range(0, no_layers):
            layer_mons = layer_distribution[l*no_layers + k]
            positions_index = sorted(set(layer_distribution)).index(layer_mons)  # The index in positions which should be read from
            for j in range(0, layer_mons):  # iterates over each monomer to create coordinates
                for i in range(0, no_ions):
                    x_values.append(b*p2p + positions[positions_index][j][no_atoms - (i + 1)][0])
                    y_values.append(c*p2p + positions[positions_index][j][no_atoms - (i + 1)][1])
                    #z_values.append(k*dist + positions[positions_index][j][no_atoms - (i + 1)][2])
                    z_values.append(k*dist + (dist/float(layer_mons))*j + positions[positions_index][j][no_atoms - (i + 1)][2])  # helix

    count = 0
    for l in range(0, no_pores):
        for k in range(0, no_layers):
            layer_mons = layer_distribution[l*no_layers + k]
            positions_index = sorted(set(layer_distribution)).index(layer_mons)  # The index in positions which should be read from
            for j in range(0, layer_mons):  # iterates over each monomer to create coordinates
                for i in range(0, no_ions):
                    count += 1
                    monomer_count += 1  # calling each sodium ion a 'monomer' for ease of numbering
                    if atom_count < 100000:
                        hundreds = int(math.floor(atom_count/100000))
                        f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(monomer_count,
                            identity[no_atoms - (no_ions - i)], identity[no_atoms - (no_ions - i)],
                            atom_count - hundreds*100000, x_values[atom_count - 1]/10.0, y_values[atom_count - 1]/10.0,
                            z_values[atom_count - 1]/10.0) + "\n")
                    atom_count += 1

    f.write('   0.00000   0.00000  0.00000\n')
    f.close()


def write_gro(positions, identity, no_layers, layer_distribution, dist, no_pores, p2p, no_ions):

    f = open('%s' % args.out, 'w')

    f.write('This is a .gro file\n')
    f.write('%s\n' % sys_atoms)

    rot = math.pi / 4

    # main monomer
    atom_count = 1
    monomer_count = 0
    for l in range(0, no_pores):  # loop to create multiple pores
        theta = 30  # angle which will be used to do hexagonal packing
        if l == 0:  # unmodified coordinates
            b = 0
            c = 0
        elif l == 1:  # move a pore directly down
            b = -1
            c = 0
        elif l == 2:  # moves pore up and to the right
            b = -math.sin(math.radians(theta))
            c = -math.cos(math.radians(theta))
        elif l == 3:  # moves a pore down and to the right
            b = math.cos(math.radians(90 - theta))
            c = -math.sin(math.radians(90 - theta))
        for k in range(no_layers):
            layer_mons = layer_distribution[l*no_layers + k]
            for j in range(layer_mons):  # iterates over each monomer to create coordinates
                monomer_count += 1
                theta = j * math.pi / (layer_mons / 2.0) + rot
                if args.offset:
                    theta += (k % 2) * (math.pi / layer_mons)
                Rx = rotate(theta)
                xyz = np.zeros(positions.shape)
                for i in range(no_atoms - no_ions):
                    if args.helix:
                        xyz[:, i] = np.dot(Rx, positions[:, i]) + [b*p2p, c*p2p, k*dist + (dist/float(layer_mons))*j]
                        hundreds = int(math.floor(atom_count/100000))
                    else:
                        xyz[:, i] = np.dot(Rx, positions[:, i]) + [b*p2p, c*p2p, k*dist]
                        hundreds = int(math.floor(atom_count/100000))
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(monomer_count, name, identity[i],
                        atom_count - hundreds*100000, xyz[0, i]/10.0, xyz[1, i]/10.0, xyz[2, i]/10.0) + "\n")
                    atom_count += 1

    # Ions:

    for l in range(no_pores):  # loop to create multiple pores
        theta = 30  # angle which will be used to do hexagonal packing
        if l == 0:  # unmodified coordinates
            b = 0
            c = 0
        elif l == 1:  # move a pore directly down
            b = - 1
            c = 0
        elif l == 2:  # moves pore up and to the right
            b = math.cos(math.radians(90 - theta))
            c = -math.sin(math.radians(90 - theta))
        elif l == 3:  # moves a pore down and to the right
            b = -math.sin(math.radians(theta))
            c = -math.cos(math.radians(theta))
        for k in range(no_layers):
            layer_mons = layer_distribution[l*no_layers + k]
            for j in range(layer_mons):  # iterates over each monomer to create coordinates
                theta = j * math.pi / (layer_mons / 2.0)
                if args.offset:
                    theta += (k % 2) * (math.pi / layer_mons)
                Rx = rotate(theta)
                xyz = np.zeros([3, no_ions])
                for i in range(0, no_ions):
                    monomer_count += 1
                    if args.helix:
                        xyz[:, i] = np.dot(Rx, positions[:, no_atoms - (i + 1)]) + [b*p2p, c*p2p, k*dist + (dist/float(layer_mons))*j]
                        hundreds = int(math.floor(atom_count/100000))
                    else:
                        xyz[:, i] = np.dot(Rx, positions[:, no_atoms - (i + 1)]) + [b*p2p, c*p2p, k*dist]
                        hundreds = int(math.floor(atom_count/100000))
                    f.write('{:5d}{:5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}'.format(monomer_count, identity[no_atoms - (i + 1)],
                        identity[no_atoms - (i + 1)], atom_count - hundreds*100000, xyz[0, i]/10.0, xyz[1, i]/10.0,
                        xyz[2, i]/10.0) + "\n")
                    atom_count += 1

    f.write('   0.00000   0.00000  0.00000\n')
    f.close()

if __name__ == "__main__":

    args = initialize()

    exec "build_mon = LC_class.%s.build_mon" % args.build_mon
    exec "no_ions = LC_class.%s.valence" % args.build_mon
    exec "name = LC_class.%s.name" % args.build_mon
    exec "plane_atoms = LC_class.%s.planeatoms" % args.build_mon
    exec "ref_index = LC_class.%s.ref_atom_index" % args.build_mon
    exec "lineatoms = LC_class.%s.lineatoms" % args.build_mon
     
    no_monomers = args.monomers  # number of monomers packed per layer around a pore
    pore_radius = args.radius  # Radius of pore (unsure of units right now)
    no_pores = args.nopores  # number of pores to be simulated
    p2p = args.p2p  # distance between pores (units tbd)
    no_layers = args.layers  # Number of layers in a pore
    dist = args.dbwl  # distance between layers (units tbd)

    f = open("%s/../Structure-Files/HII_Monomer_Configurations/%s" % (location, build_mon), "r")
    if build_mon.endswith('.pdb'):
        xyz, identity, no_atoms, lines_of_text = read_pdb_coords(f)
    elif build_mon.endswith('.gro'):
        xyz, identity, no_atoms, lines_of_text = read_gro_coords(f)
    else:
        print 'Please input a valid file type (.gro or .pdb)'

    layer_distribution = layer_dist(args.layers, args.nopores, args.layer_distribution)

    sys_atoms = sum(layer_distribution)*no_atoms  # total number of atoms in the system

    # Now rotate plane to align with xy plane
    R = rotateplane(plane_atoms, angle=(args.tilt*math.pi / 180))

    b = np.ones([1])
    for i in range(np.shape(xyz)[1]):
        coord = np.concatenate((xyz[:, i], b))
        x = np.dot(R, coord)
        xyz[:, i] = x[:3]

    # Now translate the structure to the origin

    translation = np.matrix([[1, 0, 0,-xyz[0, ref_index]], [0, 1, 0,-xyz[1, ref_index]],
                             [0, 0, 1, -xyz[2, ref_index]], [0, 0, 0, 1]])

    b = np.ones([1])
    for i in range(np.shape(xyz)[1]):
        coord = np.concatenate((xyz[:, i], b))
        x = np.dot(translation, coord)
        xyz[:, i] = x[0, :3]


    # Now rotate the xy coordinates so that the molecule is pointing towards the origin

    pt1 = [xyz[0, lineatoms[0]], xyz[1, lineatoms[0]]]  # location of C
    pt2 = [xyz[0, lineatoms[1]], xyz[1, lineatoms[1]]]  # location of C3

    origin = [0, 0]

    # find slope between two points

    m1 = slope(pt1, pt2)

    m2 = 0  # slope of line y = 0

    # find angle between lines

    theta = -math.atan((m1 - m2)/(1 + m1*m2))

    vx, vy = transdir(pt1)

    # Translation matrix
    translation = np.matrix([[1, 0, 0, vx*pore_radius*math.cos(theta)], [0, 1, 0, vy*pore_radius*math.sin(theta)],\
                             [0, 0, 1, 0], [0, 0, 0, 1]])

    b = np.ones([1])
    for i in range(np.shape(xyz)[1]):
        coord = np.concatenate((xyz[:, i], b))
        x = np.dot(translation, coord)
        xyz[:, i] = x[0, :3]

    positions = []

    for i in range(len(np.unique(layer_distribution))):  # add a list to positions for each unique value of monomers per layer
        positions.append([])
    for i in range(len(positions)):
        for j in range(0, np.sort(np.unique(layer_distribution))[i]):
            positions[i].append([])

    x_values = []  # will hold x values in the order that they appear in the positions matrix
    y_values = []  # will hold y values in the order that they appear in the positions matrix
    z_values = []  # will hold z values in the order that they appear in the positions matrix

    # rotate coordinates and store each rotated coordinate as a separate list:

    for j in range(0, len(positions)):
        for k in range(np.shape(xyz)[1]):
            x = np.array(xyz[:, k])
            for i in range(1, len(positions[j]) + 1):
                theta = i * math.pi / (len(positions[j]) / 2.0)  # angle to rotate about axis determined from no of monomers per layer
                Rx = rotate(theta)
                rot = np.dot(Rx, x)  # multiplies atomic coordinates by the rotation vector to generate new coordinates
                rot = [float(rot[0]), float(rot[1]), float(rot[2])]  # converts matrix to a list of floats
                positions[j][i - 1].append(rot)  # appends the atomic coordinates to 'positions'

    #write_gro_bak(positions, identity, no_layers, layer_distribution, dist, no_pores, p2p, no_ions)
    write_gro(xyz, identity, no_layers, layer_distribution, dist, no_pores, p2p, no_ions)
