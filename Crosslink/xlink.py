#!/usr/bin/env python
# Script to crosslink LLC monomers based on distance between carbons at the end of a simulation

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from past.utils import old_div
import os
import argparse
import numpy as np
import math
import time
import genpairs
import copy
import mdtraj as md
from llclib import transform
from matplotlib import path
import LC_class
import random

start = time.time()  # For informational purposes
location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in


def initialize():

    parser = argparse.ArgumentParser(description = 'Crosslink LLC structure')  # allow input from user

    parser.add_argument('-i', '--input', default='wiggle_init.gro', help = 'Name of input file')
    parser.add_argument('-b', '--build_mon', default='NAcarb11Vd', type=str, help='Type of monomer the system is built with')
    parser.add_argument('-l', '--layers', default=20, help='Number of layers')
    parser.add_argument('-p', '--pores', default=4, help='Number of pores')
    parser.add_argument('-n', '--no_monomers', default=6, help='Number of monomers per layer')
    parser.add_argument('-c', '--cutoff', default=5, help='Cutoff distance for cross-linking. Bottom x % of the distribution'
                                                          'of distances will be cross-linked ')
    parser.add_argument('-e', '--term_prob', default=5, type=float, help='Termination probability (%)')
    parser.add_argument('-y', '--topology', default='crosslinked_new.itp', help='Topology file that will be analyzed and modified')
    parser.add_argument('-r', '--iteration', default=0, type=int, help='Iteration number of crosslinking process')
    parser.add_argument('-d', '--cutoff_rad', default=10, help='Cutoff Distance for radical reaction')
    parser.add_argument('-x', '--xlinks', default=0, help='Total number of c1-c2 crosslinks')
    parser.add_argument('-S', '--stop', default='no', help='Is crosslinking reaction finished')
    parser.add_argument('-m', '--monomer', default='monomer2', help='Which monomer was the structure built with')
    parser.add_argument('--nogap', help='System is prepared with no gap between membrane layers', action="store_true")

    args = parser.parse_args()

    return args


def exclude_adjacent(mat_dim, images):

    # Create a matrix of 1's and 0's. Entries that are 1's should not be counted in the distance calculations because they
    # represent distances that are either between carbons on the same monomer, or are terminated carbons

    atoms = old_div(mat_dim,images)
    ones = np.ones((1, mat_dim))[0]
    exclude = np.diag(ones, 0)

    for i in range(1, images):
        ones = np.ones((1, mat_dim - i*atoms))[0]
        diag = np.diag(ones, -i*atoms)
        exclude += diag
        exclude += np.diag(ones, i*atoms)

    return exclude


def restrict(box, atom_list, buffer):

    # Get box vertices. Order is important!!! It needs to make a path that when traced from vertex to vertex is a box
    v1 = [box[0, 0, 0], box[0, 0, 1]]
    v2 = [box[0, 0, 0] + box[0, 1, 0], box[0, 1, 1]]
    v3 = [box[0, 1, 0], box[0, 1, 1]]
    v4 = [0, 0]

    # expand box isotropically according to buffer amount
    xb = box[0, 0, 0] * buffer  # x buffer
    yb = box[0, 1, 1] * buffer  # y buffer
    expand_box = np.matrix([[xb, -yb], [xb, yb], [-xb, yb], [-xb, -yb]])
    xyverts = np.array([v1, v2, v3, v4]) + expand_box

    # Define z limits
    zmax = max(box[0, 2, :])
    zmin = min(box[0, 2, :])
    thick = zmax - zmin
    zmax += thick * buffer
    zmin -= thick * buffer

    # create matplotlib path object that defines the unit cell in 2D
    p = path.Path(xyverts)

    restricted = []  # atoms that are within the distance search region

    # check all the locations to see if they are within the expanded box
    for i in range(atom_list.shape[1]):
        pt = [atom_list[0, i], atom_list[1, i]]
        if p.contains_point(pt) and zmax > atom_list[2, i] > zmin:
            restricted.append(i)

    return np.array(restricted)


def exclude(restrictions, exclusion_indices, images, natoms):
    """
    :param restrictions: existing restrictions (created using restrict function)
    :param exclusion_indices: indices of carbons that need to be excluded
    :param images: the number of periodic images. (with 3D periodic boundaries, this is 27)
    :param natoms: the total number of atoms in the unit cell
    :return: a new list of restricted atoms that does not include the excluded atoms. This list is meant to be
     passed into the distance search (genpairs.calc_dist)
    """

    nex = len(exclusion_indices)  # the number of indices to be excluded
    ex = np.zeros([nex*images])

    for j in range(images):
        count = 0
        for i in exclusion_indices:
            ex[count + j*nex] = exclusion_indices[count] + j*natoms
            count += 1

    new_restrictions = []
    for i in restrictions:
        if i not in ex:
            new_restrictions.append(i)

    return np.array(new_restrictions)


def rad_index(c1_rad_ndx_prev, c2_rad_ndx_prev, c1, c2):  # function to modify c1_rad_ndx and c2_rad_ndx since it might
    # change during manipulations of c1 and c2 lists
    c1_rad_ndx = []
    c2_rad_ndx = []
    for i in range(0, len(c1_rad_ndx_prev)):
        if c1_rad_ndx_prev[i] in c1 and c2_rad_ndx_prev[i] in c2:
            if c1_rad_ndx_prev[i] not in c1_rad_ndx and c2_rad_ndx_prev[i] not in c2_rad_ndx:
                c1_rad_ndx.append(c1_rad_ndx_prev[i])
                c2_rad_ndx.append(c2_rad_ndx_prev[i])
    return c1_rad_ndx, c2_rad_ndx


def convert_index(c1_atoms, c2_atoms, c1, c2):

    convert = []

    c1 %= 1440
    c2 %= 1440

    mon_c1 = int(old_div(c1, 3))
    mon_c2 = int(old_div(c2, 3))

    r_c1 = int(((old_div(c1, 3)) - mon_c1)*3)
    r_c2 = int(((old_div(c2, 3)) - mon_c2)*3)

    return [mon_c1*143 + int(c1_atoms[r_c1]), mon_c2*143 + int(c2_atoms[r_c2])]


def uniq(input1, input2):
    output1 = []
    output2 = []
    for i in range(0, len(input1)):
        if input1[i] not in output1 and input2[i] not in output2:
            output1.append(input1[i])
            output2.append(input2[i])
    return output1, output2

if __name__ == "__main__":

    args = initialize()

    # Read in coordinate file from simulation

    # f = open(args.input, 'r')  #open the file which has just finished being simulated and copy in all of the information
    # a = []
    # lines = 0
    # for line in f:
    #     a.append(line)
    #     lines += 1  # count number of lines in the file

    t = md.load(args.input)  # load in .gro file using mdtraj
    pos = t.xyz  # the positions from the .gro file
    box = t.unitcell_vectors

    # get all of this information from the class defining the specified build monomer
    exec("atoms = LC_class.%s.atoms" % args.build_mon)
    exec("c1_atoms = LC_class.%s.c1_atoms" % args.build_mon)
    exec("c2_atoms = LC_class.%s.c2_atoms" % args.build_mon)
    topology = "%s.itp" % args.build_mon
    exec("xlink_atoms = LC_class.%s.no_vsites" % args.build_mon)
    exec("no_dummies = LC_class.%s.no_vsites" % args.build_mon)
    exec("images = LC_class.%s.images" % args.build_mon)  # total periodic images used for distance calculations
    exec("tails = LC_class.%s.tails" % args.build_mon)

    if args.nogap:
        images *= 3

    tot_atoms = atoms*args.layers*args.no_monomers*args.pores
    tot_monomers = args.layers*args.no_monomers*args.pores


    # Coordinates of carbon 1 (the carbon on the end of the chain) and 2 (the second carbon from the end of the chains)
    C1 = np.zeros([tot_monomers*tails, 3])
    C2 = np.zeros(C1.shape)
    C1_indices = [a.index for a in t.topology.atoms if a.name in c1_atoms]  # get indices of all c1 atoms
    C2_indices = [a.index for a in t.topology.atoms if a.name in c2_atoms]  # get indices of all c2 atoms
    C1_mon_indices = np.array(C1_indices[:tails]) + 1  # add one because GROMACS counts starting at 1, not zero
    C2_mon_indices = np.array(C2_indices[:tails]) + 1

    C1t = t.atom_slice(C1_indices)  # create a new mdtraj trajectory object including only c1 atoms
    C2t = t.atom_slice(C2_indices)  # create a new mdtraj trajectory object including only c2 atoms
    C1 = C1t.xyz  # get just the positions
    C2 = C2t.xyz  # get just the positions

    no_carbons = C1.shape[1]  # number of carbon atoms that may be involved in cross-linking
    C1_grid = transform.pbcs(C1, 1, 60, box, 0, nogap=args.nogap)
    C2_grid = transform.pbcs(C2, 1, 60, box, 0, nogap=args.nogap)
    C1 = np.reshape(C1_grid, (3, C1_grid.shape[1]*C1_grid.shape[2]))
    C2 = np.reshape(C2_grid, (3, C2_grid.shape[1]*C2_grid.shape[2]))

    # Instead of calculating all distances, it is much faster to calculate only distances that are close to the boundary
    C1_restricted = restrict(box, C1, 0.1)
    C2_restricted = restrict(box, C2, 0.1)

    indices = convert_index(C1_mon_indices, C2_mon_indices, C1_restricted[0], C2_restricted[0])

    stop1 = time.time()
    print('PBCs set up: %s seconds' % (stop1 - start))
    # Find distance between carbon 1 and carbon 2 for all pairs

    mat_dim = C1.shape[1]  # Dimensions of the matrix created in the next step
    # dist = np.zeros((len(C2x), len(C1x)))

    # We need to add additional exclusions if there are any terminated/reacted atoms which should be excluded

    if int(args.iteration) != 0:

        with open(args.topology, 'r') as f:  # take a look at the topology from the previous iteration

            top = []  # put all of the lines into a list
            for line in f:
                top.append(line)

        # look through each line for the 'T' marker and an '*' marker
        exclude_T_c1 = []
        exclude_T_c2 = []
        reactive_c2 = []
        c1_ndx = []  # we also need to number c1 and c2 in the context of other c1's and c2's
        c2_ndx = []  # i.e. we need numbers of c1-1, c1-2 ... and c2-1, c2-2 ... and see how they relate to the indices

        # find the [ atoms ] section
        atoms_index = 0  # find index where [ atoms ] section begins
        while top[atoms_index].count('[ atoms ]') == 0:
            atoms_index += 1

        atoms_end = atoms_index + 2  # start looking on a line where there isn't text
        while top[atoms_end] != '\n':  # look through the whole [ atoms ] section
            if str.strip(top[atoms_end][22:28]) in c1_atoms:  # make a list of the indices of all c1 atoms
                c1_ndx.append(int(top[atoms_end][0:5]))
            if str.strip(top[atoms_end][22:28]) in c2_atoms:  # make a list of the indices of all c2 atoms
                c2_ndx.append(int(top[atoms_end][0:5]))
            if str.strip(top[atoms_end])[-1] == 'T':  # see which lines are terminated and should therefore be excluded
                if str.strip(top[atoms_end][22:28]) in c1_atoms:  # if its a c1, add it to a list of terminated c1's
                    exclude_T_c1.append(int(top[atoms_end][0:5]))
                if str.strip(top[atoms_end][22:28]) in c2_atoms:  # if its a c2, add it to a list of terminated c2's
                    exclude_T_c2.append(int(top[atoms_end][0:5]))
            if str.strip(top[atoms_end])[-1] == '*':  # see which atoms are reactive radicals
                reactive_c2.append(int(top[atoms_end][0:5]))  # if it is, add it too a list of reactive c2 atoms
            atoms_end += 1  # increment while loop

        c1_ex_no = []  # Now number the c1's out of all c1's
        c2_ex_no = []  # Same as above but for c2
        c2_reactive_no = []  # and the radical
        for i in range(len(exclude_T_c1)):
            c1_ex_no.append(c1_ndx.index(exclude_T_c1[i]))  # their number is based on the index where they are located
        for i in range(len(exclude_T_c2)):  #exclude_T_c2 and exclude_T_c1 are not necessarily the same length so we need a second loop
            c2_ex_no.append(c2_ndx.index(exclude_T_c2[i]))
        for i in range(len(reactive_c2)):
            c2_reactive_no.append(c2_ndx.index(reactive_c2[i]))  # the indices of the radical carbons

        # for i in c1_ex_no:  # now fill the exclusions matrix with 1's where we don't want a distance calc.
        #     for j in range(0, images):  # do it for all of the periodic cells
        #         exclude[:, j*mat_dim/images + i] = 1  # excludes across all pbc's too
        #
        # for i in c2_ex_no:  # do the same for c2
        #     for j in range(0, images):
        #         exclude[j*mat_dim/images + i, :] = 1

        C1_restricted = exclude(C1_restricted, c1_ex_no, images, old_div(C1.shape[1],images))
        C2_restricted = exclude(C2_restricted, c2_ex_no, images, old_div(C2.shape[1],images))

    print('beginning distance calculation')

    start_dist = time.time()
    dist = genpairs.calc_dist2(C1, C2, C1_restricted, C2_restricted)

    stop3 = time.time()

    d = [d for d in dist[:, C2_restricted[1]] if d < 1000]

    indices = convert_index(C1_mon_indices, C2_mon_indices, C1_restricted[np.argmin(d)], C2_restricted[1])

    print('Distances Calculated: %s seconds' % (stop3 - start_dist))

    # Now see which of these distances meet the cutoff criteria
    # Find the distance of the closest carbon to each c1 and its index

    ncarb1 = C1.shape[1]
    ncarb2 = C2.shape[1]
    min_dist = np.zeros([ncarb1]) + 1000
    min_index = np.zeros([ncarb1])  # index of minimum value of distances for each monomer-monomer measurement
    for i in C1_restricted:
        min_dist[i] = min(dist[i, :])
        min_index[i] = np.argmin(dist[i, :]) % (old_div(ncarb1,images))  # index corresponds to the monomer with which

    change_index1 = []  # index of C1 from min_dist which needs to be changed
    change_index2 = []  # index of C2 from dist which needs to be changed

    if int(args.iteration) != 0:

        # Looking just at the radical reactive sites. Does the same thing as the previous block but for a smaller list

        c2_radicals = []
        for i in C2_restricted:
            if i % (old_div(ncarb1,images)) in c2_reactive_no:
                c2_radicals.append(i)

        min_dist_rad = np.zeros([len(c2_radicals)])
        min_index_rad = np.zeros([len(c2_radicals)])
        for i in range(len(c2_radicals)):
            min_dist_rad[i] = min(dist[c2_radicals[i], :])
            min_index_rad[i] = np.argmin(dist[c2_radicals[i], :]) % (old_div(ncarb1,images))

        min_list_rad = min_dist_rad.tolist()
        for i in range(int((old_div(float(args.cutoff_rad),100))*len(min_list_rad))):  # This could be done by sorting then cutting off first ten percent -- its probably better that way
            m = min(min_list_rad)
            min_list_rad.remove(m)

        if min_list_rad:
            cutoff = min(min_list_rad)  # Now the cutoff value is the minimum of min_list_rad
        else:
            cutoff = 0

        if cutoff > 0.7:
            cutoff = 0.7

        for i in range(len(min_index_rad)):  # Add to change_index if they meet the cutoff criteria
            if min_dist_rad[i] < cutoff:
                change_index1.append(int(min_index_rad[i]))
                change_index2.append(c2_radicals[i])

        no_radical_rxns = len(change_index2)

    # find the cutoff distance for cross-linking

    min_list = [a for a in min_dist if a < 1000]

    distances = []
    for i in range(int((old_div(float(args.cutoff),100))*len(min_list))):  # looks at a percentage of the total values based on user input of cutoff
        m = min(min_list)  # finds minimum of min_list
        distances.append(m)
        min_list.remove(m)  # removes that value from min_list

    cutoff = min(min_list)  # The minimum value left after the modification of min_list is the cutoff value
    print('Cutoff = %s nm' % cutoff)
    if cutoff > 0.6:
        cutoff = 0.6

    count = 0
    for i in C1_restricted:
        if min_dist[i] < cutoff:  # find distances which meet the cutoff criteria
            change_index1.append(i % (old_div(ncarb1,images)))  # holds the index of the atom associated with the met criteria for primary carbons (C20, C34, C48)
            change_index2.append(int(min_index[i]))  # Same as above but for the secondary carbons
            count += 1

    # Now that everything has an index that needs to be changed, we must interpret those indices

    # C19 is the 26th atom, C20 is 27th, C33 is 41st, C34 is 42nd, C47 is 56th, C48 is 57th, in each monomer

    # We also need to know which index refers to which monomer and tail -- applies equally for C1 and C2
    # Every third index starts a new monomer. Each index in between is a tail (inclusive)

    tail1 = np.arange(0, ncarb1 + 1, 3)  # indices of carbons
    tail2 = np.arange(1, ncarb1 + 1, 3)
    tail3 = np.arange(2, ncarb1 + 1, 3)

    C1_no = []
    C2_no = []

    for i in range(len(change_index1)):
        if change_index1[i] in tail1:  # This is C1 therefore if this is true, then the atom is C20 (atom no 27)
            C1_no.append((old_div(change_index1[i],3))*atoms + C1_mon_indices[0])
        if change_index1[i] in tail2:  # This is C1 therefore if this is true, then the atom is C34 (atom no 42)
            C1_no.append((old_div(change_index1[i],3))*atoms + C1_mon_indices[1])  # division automatically round down
        if change_index1[i] in tail3:  # This is C1 therefore if this is true, then the atom is C48 (atom no 57)
            C1_no.append((old_div(change_index1[i],3))*atoms + C1_mon_indices[2])  # division automatically round down
        if change_index2[i] in tail1:  # This is C2 therefore if this is true, then the atom is C19 (atom no 26)
            C2_no.append((old_div(change_index2[i],3))*atoms + C2_mon_indices[0])
        if change_index2[i] in tail2:  # This is C2 therefore if this is true, then the atom is C33 (atom no 41)
            C2_no.append((old_div(change_index2[i],3))*atoms + C2_mon_indices[1])  # division automatically round down
        if change_index2[i] in tail3:  # This is C2 therefore if this is true, then the atom is C47 (atom no 56)
            C2_no.append((old_div(change_index2[i],3))*atoms + C2_mon_indices[2])  # division automatically round down

    if int(args.iteration) != 0:

        c1_rad_ndx = []
        c2_rad_ndx = []
        for i in range(no_radical_rxns):  # The first entries in C2_no are all radicals up to no_radical_rxns
            c2_rad_ndx.append(C2_no[i])
            c1_rad_ndx.append(C1_no[i])

        c1_rad_ndx_prev = c1_rad_ndx
        c2_rad_ndx_prev = c2_rad_ndx

    c1_uniq, c2_uniq = uniq(C1_no, C2_no)

    if int(args.iteration) != 0:
        c1_rad_ndx, c2_rad_ndx = rad_index(c1_rad_ndx_prev, c2_rad_ndx_prev, c1_uniq, c2_uniq)

    # In some cases, it is possible that a c1 and a c2 on the same chain may be involved in independent cross-link reactions
    # We need to get rid of those to avoid unnecessary confusion. Since c2 is the reactive atom (it is a radical at some
    # point), it will trump c1. So if a c2 and c1 are on the same chain, then c1 will be removed from consideration.

    # Create a list of the indices of all the c1 atoms adjacent to reacting c2's

    if int(args.iteration) != 0:  # preserve previous radical indices for testing later since some may be eliminated
        c1_rad_ndx_prev = c1_rad_ndx
        c2_rad_ndx_prev = c2_rad_ndx

    c1 = []
    c2 = []

    for i in range(len(c1_uniq)):
        if c2_uniq[i] + 1 in c1_uniq:
            if c2_uniq[i] + 1 not in c1:
                c1.append(c1_uniq[i])
                c2.append(c2_uniq[i])
        else:
            c1.append(c1_uniq[i])
            c2.append(c2_uniq[i])

    # print("c1: %s" %c1)
    # print("c2: %s" %c2)
    # A table of crosslinks
    print('Crosslinks:')
    print('{:^10s}{:^10s}'.format('c1', 'c2'))
    print('--------------------')
    for i in range(len(c1)):
        print('{:^8d} -- {:^8d}'.format(c1[i], c2[i]))


    # adjacent_c1 = []
    # for i in range(len(c2_uniq)):
    #     adjacent_c1.append(c2_uniq[i] + 1)
    #
    # # Create the final lists of c1 and c2 that will be cross-linked
    # c1 = []
    # c2 = []
    # for i in range(len(adjacent_c1)):
    #     if c1_uniq[i] not in adjacent_c1:  # add values to the lists that have c1's that are not adjacent to reacting c2's
    #         c1.append(c1_uniq[i])
    #         c2.append(c2_uniq[i])

    if int(args.iteration) != 0:
        c1_rad_ndx, c2_rad_ndx = rad_index(c1_rad_ndx_prev, c2_rad_ndx_prev, c1, c2)

    # Some reactions will terminate instead of cross linking. This will cause the dummy hydrogen on C2 to become active

    if int(args.iteration) != 0:  # preserve our old list of radicals again
        c1_rad_ndx_prev = c1_rad_ndx
        c2_rad_ndx_prev = c2_rad_ndx

    tp = args.term_prob  # termination probability

    # indices = random.sample(list(range(0, int(100/tp))), 1)
    term_prob_array = np.zeros([int(100/tp)])  # make an array with a a percentage of 1's equal to the termination probability
    # for i in indices:
    #     term_prob_array[i] = 1
    term_prob_array[random.randint(0, int(100/tp) - 1)] = 1

    # Termination occurs at C2
    # Move termination carbons to a new list
    term_c1 = []
    term_c2 = []
    term_indices = []
    no_xlinks = copy.deepcopy(len(c2))  # The length of c2 may change in the next loops so I am preserving it here

    # Add the termination carbons to new list
    for i in range(no_xlinks):
        term = np.random.choice(term_prob_array)
        if term == 1:
            term_c1.append(c1[i])
            term_c2.append(c2[i])
            term_indices.append(i)

    # Remove terminated carbons from c1 and c2

    for i in range(len(term_indices)):
        len_diff = no_xlinks - len(c1)  # Account for varying size of c1 and c2 lists
        del c1[term_indices[i] - len_diff]
        del c2[term_indices[i] - len_diff]

    if int(args.iteration) != 0:
        c1_rad_ndx, c2_rad_ndx = rad_index(c1_rad_ndx_prev, c2_rad_ndx_prev, c1, c2)  # rewrite c1_rad_ndx and c2_rad_ndx
        # Separate the termination lists into a list for termination immediately upon initiation (term) and a list for
        # termination occurring on already initiated radicals
        reactive_c2_term = []  # reactive c2's which didn't meet distance cut-off criteria but, since there is still a
        # possibility for these to terminate, we must account for them
        for i in reactive_c2:  # This list contains all the entries in c2_rad_ndx
            if i not in c2_rad_ndx:  # so first make sure that we aren't adding something from c2_rad_ndx to the list
                term = np.random.choice(term_prob_array)
                if term == 1:
                    reactive_c2_term.append(i)  # These carbons will terminate
        term_rad_c1 = []  # radicals from previous iterations which will terminate
        term_rad_c2 = []
        for i in term_c2:  # look at all the entries in term
            if i in c2_rad_ndx:  # if it is in c2_rad_ndx
                term_rad_c1.append(term_c1.index(i))  # then find the corresponding c1 and record it
                term_rad_c2.append(i)  # record that value (c1 and c2 are in pairs)
                del term_c1[term_c1.index(i)]  # now delete those values from the previous term lists
                del term_c2[term_c2.index(i)]

        if args.stop == 'yes':
            # find the [ atoms ] section ... again
            atoms_index = 0  # find index where [ atoms ] section begins
            while top[atoms_index].count('[ atoms ]') == 0:
                atoms_index += 1

            atoms_end = atoms_index + 2  # start looking on a line where there isn't text
            count = 0
            while top[atoms_end] != '\n':  # look through the whole [ atoms ] section
                # print str.strip(top[atoms_end])[-1]
                if str.strip(top[atoms_end])[-1] == '*':  # see which atoms are reactive radicals
                    count += 1
                    if int(top[atoms_end][0:5]) not in reactive_c2_term:
                        reactive_c2_term.append(int(top[atoms_end][0:5]))  # if it is, add it too a list of reactive c2 atoms
                atoms_end += 1

    # For the first iteration, use Assembly_itp.py to create a topology for the entire assembly and then write it to a file
    # which will be edited to incorporate cross-links. Otherwise, read the topology from the previous iteration

    if int(args.iteration) == 0:
        # import subprocess
        # import os
        # location = os.environ['GITHUB']  # if there is an error here, you need to add the path to where all of the github
        # # files are stored to an environment variable in your .bashrc
        # subprocess.call(["python", "%s/Scripts/Assembly_itp.py" % location, "-x", "on", "-m", "%s" % args.monomer,
        #                  "-O", "crosslinked.itp"])
        #
        # # open and read that new file
        from llclib import file_rw

        file_rw.write_assembly(args.build_mon, 'on', "crosslinked.itp", ncarb1/tails/images)
        print('Crosslinked.itp file written')

        f = open('crosslinked.itp', 'r')
        b = []
        for line in f:
            b.append(line)
    else:
        f = open(args.topology, 'r')
        b = []
        for line in f:
            b.append(line)
        print('%s read into list' % args.topology)

    # Make lists with numbers of H atoms whose type needs to change as a consequence of the cross linking reaction

    # # first, redefine tail 1, 2 and 3 using C2 atom numbers as a reference
    # # find the residue numbers of the first C33, C47, C19 (C2 carbons)
    #
    # index = 0
    # #while str.strip(b[index][22:29]) != 'C19':
    # while str.strip(b[index][23:28] not in c1_atoms == 0:
    #     index += 1  # increments the while loop
    #
    # C2_1 = int(b[index][0:5])
    #
    # index = 0
    # while b[index].count('C33') == 0:
    # #while str.strip(b[index][22:29]) != 'C33':
    #     index += 1  # increments the while loop
    #
    # C2_2 = int(b[index][0:5])
    # index = 0
    #
    # while b[index].count('C47') == 0:
    # #while str.strip(b[index][22:29]) != 'C47':
    #     index += 1  # increments the while loop
    #
    # C2_3 = int(b[index][0:5])

    index = 0
    #while str.strip(b[index][22:29]) != 'C19':
    while b[index].count('C19') == 0:
        index += 1  # increments the while loop

    C2_1 = int(b[index][0:5])

    index = 0
    while b[index].count('C33') == 0:
    #while str.strip(b[index][22:29]) != 'C33':
        index += 1  # increments the while loop

    C2_2 = int(b[index][0:5])
    index = 0

    while b[index].count('C47') == 0:
    #while str.strip(b[index][22:29]) != 'C47':
        index += 1  # increments the while loop

    C2_3 = int(b[index][0:5])

    tail1 = []
    tail2 = []
    tail3 = []

    for i in range(0, tot_monomers):
        tail1.append(atoms*i + C2_1)
        tail1.append(atoms*i + C2_1 + 1)
        tail2.append(atoms*i + C2_2)
        tail2.append(atoms*i + C2_2 + 1)
        tail3.append(atoms*i + C2_3)
        tail3.append(atoms*i + C2_3 + 1)

    # tail1 = []
    # tail2 = []
    # tail3 = []
    #
    # for i in range(0, tot_monomers):
    #     tail1.append(atoms*i + C2_mon_indices[0])  # this can obviously be improved to the commented loop below
    #     tail1.append(atoms*i + C1_mon_indices[0])
    #     tail2.append(atoms*i + C2_mon_indices[1])
    #     tail2.append(atoms*i + C1_mon_indices[1])
    #     tail3.append(atoms*i + C2_mon_indices[2])
    #     tail3.append(atoms*i + C1_mon_indices[2])

    # alltails = np.zeros([tails, 2*tot_monomers])
    # for i in range(tails):
    #     for j in range(tot_monomers):
    #         alltails[i, 2*j] = C2_mon_indices[i]
    #         alltails[i, 2*j + 1] = C1_mon_indices[i]

    # We need to change the atom type of the non-bonding c1. If it is initiated, then it goes from c2 to c3
    # We also need a list of c2's which are left as radicals. They are adjacent to c1
    other_c1 = []
    radical_c2 = []
    for i in range(0, len(c2)):
        other_c1.append(c2[i] + 1)
        radical_c2.append(c1[i] - 1)

    if int(args.iteration) != 0:
        # break apart the c1 and c2 list for now while atom types are being changed
        del c1[:len(c1_rad_ndx)]
        del c2[:len(c2_rad_ndx)]

    H = []  # list of residue numbers of H that need to be changed from ha to hc -- unneeded. Delete once this whole thing works
    H_new1 = []  # list of residue number of H's that will be converted from dummy atoms to real atoms
    for i in range(0, len(c2)):
        if c2[i] in tail1:  # i.e. C19
            H.append(c2[i] + 59)
            H.append(c2[i] + 60)
            H.append(c2[i] + 61)
            if c1[i] in tail1:
                H.append(c1[i] + 59)
                H.append(c1[i] + 60)
            if c1[i] in tail2:
                H.append(c1[i] + 69)
                H.append(c1[i] + 70)
            if c1[i] in tail3:
                H.append(c1[i] + 79)
                H.append(c1[i] + 80)
            H_new1.append(c2[i] + 113)
        if c2[i] in tail2:  # i.e. C33
            H.append(c2[i] + 69)
            H.append(c2[i] + 70)
            H.append(c2[i] + 71)
            if c1[i] in tail1:
                H.append(c1[i] + 59)
                H.append(c1[i] + 60)
            if c1[i] in tail2:
                H.append(c1[i] + 69)
                H.append(c1[i] + 70)
            if c1[i] in tail3:
                H.append(c1[i] + 79)
                H.append(c1[i] + 80)
            H_new1.append(c2[i] + 100)
        if c2[i] in tail3:  # i.e. C47
            H.append(c2[i] + 79)
            H.append(c2[i] + 80)
            H.append(c2[i] + 81)
            if c1[i] in tail1:
                H.append(c1[i] + 59)
                H.append(c1[i] + 60)
            if c1[i] in tail2:
                H.append(c1[i] + 69)
                H.append(c1[i] + 70)
            if c1[i] in tail3:
                H.append(c1[i] + 79)
                H.append(c1[i] + 80)
            H_new1.append(c2[i] + 87)

    # need to do the same for the termination carbons except we need to include the dummy hydrogen and don't have to worry
    # about the term_c1 list since nothing is happening on that end anymore

    # the c2 where termination is happening is no longer the same c2 which is referenced in term_c2
    term = []  # define those termination atoms as 'term' which are the carbons adjacent to term_c1
    for i in range(len(term_c1)):
        term.append(term_c1[i] - 1)

    H_new2 = []  # keep these separate for indexing purposes. In H_new2, termination is based on the tail containing c1
        # since initiation occurs at c1.

    for i in range(len(term)):
        if term[i] in tail1:  # i.e. C19
            H.append(term[i] + 59)
            H.append(term[i] + 60)
            H.append(term[i] + 61)
            H_new2.append(term[i] + 112)
            H_new2.append(term[i] + 113)
        if term[i] in tail2:  # i.e. C33
            H.append(term[i] + 69)
            H.append(term[i] + 70)
            H.append(term[i] + 71)
            H_new2.append(term[i] + 99)
            H_new2.append(term[i] + 100)
        if term[i] in tail3:  # i.e. C47
            H.append(term[i] + 79)
            H.append(term[i] + 80)
            H.append(term[i] + 81)
            H_new2.append(term[i] + 86)
            H_new2.append(term[i] + 87)

    if int(args.iteration) != 0:
        H_new3 = []  # termination at radicals that exist from previous iterations will occur at c2. The terminated chain is
        # determined by the location of the radical c2. Nothing happens on the tail containing c1
        H_rad_carbons = term_rad_c2 + reactive_c2_term
        for i in range(0, len(H_rad_carbons)):
            if H_rad_carbons[i] in tail1:
                H_new3.append(H_rad_carbons[i] + 112)
            if H_rad_carbons[i] in tail2:
                H_new3.append(H_rad_carbons[i] + 99)
            if H_rad_carbons[i] in tail3:
                H_new3.append(H_rad_carbons[i] + 86)

    # find the indices of all fields that need to be modified

    # [ atoms ]

    atoms_index = 0  # find index where [ atoms ] section begins
    while b[atoms_index].count('[ atoms ]') == 0:
        atoms_index += 1

    atoms_count = atoms_index + 2
    while b[atoms_count] != '\n':
        atoms_count += 1  # increments the while loop

    if args.iteration != 0:
        if c1 == [] and c2 == [] and other_c1 == [] and radical_c2 == [] and c1_rad_ndx == [] \
                    and c2_rad_ndx == [] and H_new1 == [] and cutoff == 0.6:
            Stop_next_iter = 1  # finish this iteration, then next time the system will be completely terminated
        else:
            Stop_next_iter = 0
    else:
        Stop_next_iter = 0

    if args.stop == 'yes':
        c1 = []
        c2 = []
        other_c1 = []
        radical_c2 = []
        c1_rad_ndx = []
        c2_rad_ndx = []
        H_new1 = []

    count = 0
    for i in range(atoms_index + 2, atoms_count):
        res_num = int(b[i][0:5])
        if res_num in c1:  # change bonding carbon 1 from c2 to c3 since it becomes sp3 hybridized
            b[i] = b[i][0:5] + b[i][5:10].replace('c2', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + ' ;T' + '\n'
        if res_num in other_c1:  # The c1 on the same tail as the bonding c2 is now sp3 hybridized
            b[i] = b[i][0:5] + b[i][5:10].replace('c2', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + ' ;T' + '\n'
        if res_num in c2:  # the bonding c2 becomes sp3 hybridized
            if b[i][5:10].count('ce') == 1:
                b[i] = b[i][0:5] + b[i][5:10].replace('ce', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + ' ;T' + '\n'
            if b[i][5:10].count('c2') == 1:
                b[i] = b[i][0:5] + b[i][5:10].replace('c2', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + ' ;T' + '\n'
        if res_num in term:  # the terminating c2 (where the termination actually happens) becomes sp3 hybridized
            b[i] = b[i][0:5] + b[i][5:10].replace('ce', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + ' ;T' + '\n'  # if it doesn't find the string to replace it does nothing
        if res_num in H_new1:  # Dummy atoms become real during the bonding process and have mass
            b[i] = b[i][0:5] + b[i][5:15].replace('hc_d', 'hc  ') + b[i][15:53] + b[i][53:61].replace('0.00000', '1.00800') + b[i][61:len(b[atoms_index + 2])]
        if res_num in H_new2:  # Dummy atoms which become real during the termination process and have mass
            b[i] = b[i][0:5] + b[i][5:15].replace('hc_d', 'hc  ') + b[i][15:53] + b[i][53:61].replace('0.00000', '1.00800') + b[i][61:len(b[atoms_index + 2])]
        if res_num in radical_c2:  # mark the c2 atoms that are now radicals. They remain sp2 hybridized but are reactive
            b[i] = b[i][0:5] + b[i][5:10].replace('ce', 'c2') + b[i][10:len(b[atoms_index + 2]) - 1] + ' ;*' + '\n'
        if res_num in term_c1:  # c1 in a terminating chain will be sp3
            b[i] = b[i][0:5] + b[i][5:10].replace('c2', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + ' ;T' + '\n'
        if int(args.iteration) != 0:
            if res_num in c2_rad_ndx or res_num in c1_rad_ndx or res_num in term or res_num in reactive_c2_term or res_num in term_rad_c2:
                b[i] = b[i][0:5] + b[i][5:10].replace('c2', 'c3') + b[i][10:len(b[atoms_index + 2]) - 1] + ' ;T' + '\n'
            if res_num in H_new3:
                b[i] = b[i][0:5] + b[i][5:15].replace('hc_d', 'hc  ') + b[i][15:53] + b[i][53:61].replace('0.00000', '1.00800') + b[i][61:len(b[atoms_index + 2])]

    if int(args.iteration) != 0:
        # re-assemble the full c1, c2 list
        c1 = c1_rad_ndx + c1
        c2 = c2_rad_ndx + c2


    # Before going further, we need to label all of the dummy atoms that are still dummies so that they are not written into
    # the bonds, angles, dihedrals or pairs section

    # generate a list of indices of all the dummy atoms including those who will change to real atoms
    dummy_list = []
    # for i in range(0, tot_monomers):
    #     for k in range(0, no_dummies):
    #         dummy_list.append(i*atoms + (atoms - k))  # dummies are listed at the end of the atomtypes
    for i in range(atoms_index + 2, atoms_count):
        if str.strip(b[i][5:13]) == 'hc_d':
            dummy_list.append(int(b[i][0:5]))

    leftovers = []  # poor dummies that didn't get turned into real atoms
    for i in range(0, len(dummy_list)):
        if dummy_list[i] not in H_new1 and dummy_list[i] not in H_new2:  # sees if dummies were turned to real atoms
            leftovers.append(dummy_list[i])  # if not, they become part of the leftovers

    # [ bonds ]

    bonds_index = 0  # find index where [ bonds ] section begins
    while b[bonds_index].count('[ bonds ]') == 0:
        bonds_index += 1

    nb = 0  # number of lines in the 'bonds' section
    bond_list = []
    bond_count = bonds_index + 2
    while b[bond_count] != '\n':
        a_b = int(b[bond_count][0:6])
        b_b = int(b[bond_count][6:13])
        if a_b not in leftovers and b_b not in leftovers:
            bond_list.append([a_b, b_b])
            nb += 1  # counting number of lines in 'bonds' section
        bond_count += 1  # increments while loop

    # properly format bonds list
    bonds = []
    for i in range(0, len(bond_list)):
        bonds.append('{:6d}{:7d}{:4d}'.format(bond_list[i][0], bond_list[i][1], 1))

    # eliminate the old bonds list

    del b[bonds_index + 2: bond_count]

    # Insert new bonds list

    count = bonds_index + 2
    for i in range(0, len(bonds)):
        b.insert(count, bonds[i] + '\n')
        count += 1

    # Add bonds between cross-linking atoms
    xlinks = int(args.xlinks)
    for i in range(0, len(c1)):
        b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(c1[i], c2[i], 1) + "\n")
        xlinks += 1
        nb += 1

    # Add new H bonds formed during propagation (dummy H's becoming real)
    if int(args.iteration) != 0:
        for i in range(0, len(H_new1)):
            b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(c2[i + len(c2_rad_ndx)] + 1, H_new1[i], 1) + "\n")
            nb += 1
    else:
        for i in range(0, len(H_new1)):
            b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(c2[i] + 1, H_new1[i], 1) + "\n")
            nb += 1

    # new H bonds formed during termination (2 dummy H's becoming real)
    k = 0

    for i in range(len(term)):
        b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(term[i], H_new2[k], 1) + "\n")
        k += 1
        nb += 1
        b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(term[i] + 1, H_new2[k], 1) + "\n")
        nb += 1
        k += 1

    if int(args.iteration) != 0:
        # new H bonds formed when radicals from previous iterations are terminated
        for i in range(0, len(H_new3)):
            b.insert(nb + bonds_index + 2, '{:6d}{:7d}{:4d}'.format(H_rad_carbons[i], H_new3[i], 1) + "\n")
            nb += 1

    # Now we need to generate dihedrals, angles and pairs lists. We have all the atoms we need. Let's make one list that
    # holds all of these values:

    if int(args.iteration) == 0:
        atoms_of_interest = radical_c2 + c1 + c2 + other_c1 + term + term_c1 + H + H_new1 + H_new2
    else:
        atoms_of_interest = radical_c2 + c1 + c2 + other_c1 + term + term_c1 + H + H_new1 + H_new2 + term_rad_c2 + H_new3 + reactive_c2_term

    # See genpairs.pyx to see what is going on in the following lines:

    # Generate a list of all bonds. 'list' contains all bonds. 'Condensed_list' contains bonds only relevant to 'atoms_of_interest'
    list, condensed_list = genpairs.read_lists(b, atoms_of_interest, bonds_index + 2, nb + bonds_index + 2)

    # Generate a pairs list, a dihedrals list, and an angles list that includes those formed by the newly created bonds

    pairs_list, dihedrals_list, angles_list = genpairs.gen_pairs_list(list, condensed_list)

    # Find where the [ pairs ] section is located

    pairs_index = 0  # find index where [ pairs ] section begins
    while b[pairs_index].count('[ pairs ]') == 0:
        pairs_index += 1

    npair = 0  # number of lines in the 'pairs' section
    pairs_count = pairs_index + 2  # keep track of index of a
    pairs_list_prev = []
    while b[pairs_count] != '\n':
        a_p = b[pairs_count][0:6]
        b_p = b[pairs_count][6:13]
        if a_p not in leftovers and b_p not in leftovers:
            pairs_list_prev.append([a_p, b_p])  # record pairs that already exist
            npair += 1
        pairs_count += 1

    # Eliminate duplicates in the entire pairs list
    unique_pairs = genpairs.uniq_pair(pairs_list, pairs_list_prev)

    # Create a properly formatted pairs list
    pairs = []
    for i in range(0, len(unique_pairs)):
        pairs.append('{:6d}{:7d}{:7d}'.format(int(unique_pairs[i][0]), int(unique_pairs[i][1]), 1))

    # eliminate the old pairs list

    del b[pairs_index + 2: pairs_count]

    # Insert new pairs list

    count = pairs_index + 2
    for i in range(0, len(pairs)):
        b.insert(count, pairs[i] + '\n')
        count += 1

    # [ angles ]

    angles_index = 0  # find index where [ angles ] section begins
    while b[angles_index].count('[ angles ]') == 0:
        angles_index += 1

    angle_count = angles_index + 2  # keep track of index of a
    angles_prev = []
    while b[angle_count] != '\n':
        a_a = int(b[angle_count][0:6])
        b_a = int(b[angle_count][6:13])
        c_a = int(b[angle_count][13:20])
        if a_a not in leftovers and b_a not in leftovers and c_a not in leftovers:
            angles_prev.append([a_a, b_a, c_a])
        angle_count += 1

    # Eliminate duplicates in the entire angles list
    unique_angles = genpairs.uniq_angles(angles_list, angles_prev)

    # Create a properly formatted angles list
    angles = []
    for i in range(0, len(unique_angles)):
        angles.append('{:6d}{:7d}{:7d}{:7d}'.format(unique_angles[i][0], unique_angles[i][1],
                                                            unique_angles[i][2], 1))

    # eliminate the old angles list

    del b[angles_index + 2: angle_count]

    # Insert new angles list

    count = angles_index + 2
    for i in range(0, len(angles)):
        b.insert(count, angles[i] + "\n")
        count += 1

    # [ dihedrals ] ; propers

    dihedrals_p_index = 0  # find index where [ dihedrals ] section begins (propers)
    while b[dihedrals_p_index].count('[ dihedrals ] ; propers') == 0:
        dihedrals_p_index += 1

    dihedrals_p_count = dihedrals_p_index + 2  # keep track of index of a
    dihedrals_prev = []
    while b[dihedrals_p_count] != '\n':
        a_d = int(b[dihedrals_p_count][0:6])
        b_d = int(b[dihedrals_p_count][6:13])
        c_d = int(b[dihedrals_p_count][13:20])
        d_d = int(b[dihedrals_p_count][20:27])
        if a_d not in leftovers and b_d not in leftovers and c_d not in leftovers and d_d not in leftovers:
            dihedrals_prev.append([a_d, b_d, c_d, d_d])
        dihedrals_p_count += 1

    # Eliminate duplicates in entire dihedrals list
    unique_dihedrals = genpairs.uniq_dihedral(dihedrals_list, dihedrals_prev)

    # properly format dihedrals list
    dihedrals = []
    for i in range(0, len(unique_dihedrals)):
        dihedrals.append('{:6d}{:7d}{:7d}{:7d}{:4d}'.format(unique_dihedrals[i][0], unique_dihedrals[i][1],
                                                            unique_dihedrals[i][2], unique_dihedrals[i][3], 3))

    # eliminate the old dihedrals list

    del b[dihedrals_p_index + 2: dihedrals_p_count]

    # Insert new dihedrals list

    count = dihedrals_p_index + 2
    for i in range(0, len(dihedrals)):
        b.insert(count, dihedrals[i] + "\n")
        count += 1

    # [ dihedrals ] ; impropers

    # In places where the double bond b/w c1 and c2 has changed to a single bond, the improper dihedrals need to be removed
    # generate a list of impropers of interest (ones that will potentially be eliminated)

    impropers_list = []
    for i in range(0, tot_monomers):
        impropers_list.append([atoms*i + 26, atoms*i + 27, atoms*i + 86, atoms*i + 87])
        impropers_list.append([atoms*i + 41, atoms*i + 42, atoms*i + 111, atoms*i + 112])
        impropers_list.append([atoms*i + 56, atoms*i + 57, atoms*i + 136, atoms*i + 137])
        impropers_list.append([atoms*i + 25, atoms*i + 26, atoms*i + 27, atoms*i + 85])
        impropers_list.append([atoms*i + 40, atoms*i + 41, atoms*i + 42, atoms*i + 110])
        impropers_list.append([atoms*i + 55, atoms*i + 56, atoms*i + 57, atoms*i + 135])

    imp_remove_c = term + term_c1 + other_c1 + c1 + c2  # carbons that are part of improper dihedrals that will be removed
    imp_of_interest = []
    for i in range(0, len(impropers_list)):
        if impropers_list[i][1] in imp_remove_c:
            imp_of_interest.append(impropers_list[i])

    dihedrals_imp_index = 0  # find index where [ dihedrals ] section begins (impropers)
    while b[dihedrals_imp_index].count('[ dihedrals ] ; impropers') == 0:
        dihedrals_imp_index += 1

    dihedrals_imp_count = dihedrals_imp_index + 2
    start_imp = dihedrals_imp_count

    dihedrals_imp = []
    while b[dihedrals_imp_count] != '\n':
        a_imp = [int(b[dihedrals_imp_count][0:6]), int(b[dihedrals_imp_count][6:13]), int(b[dihedrals_imp_count][13:20]), int(b[dihedrals_imp_count][20:27])]
        a_imp.sort()
        if a_imp not in imp_of_interest:
            dihedrals_imp.append(b[dihedrals_imp_count])
        dihedrals_imp_count += 1

    # eliminate the old improper dihedrals list

    del b[dihedrals_imp_index + 2: dihedrals_imp_count]

    # Insert new improper dihedrals list

    count = dihedrals_imp_index + 2
    for i in range(0, len(dihedrals_imp)):
        b.insert(count, dihedrals_imp[i])
        count += 1

    # [ virtual_sites4 ]

    vsite_index = 0  # find index where [ virtual_sites4 ] section begins (propers)
    while b[vsite_index].count('[ virtual_sites4 ]') == 0:
        vsite_index += 1

    # vsite_count = vsite_index + 2

    # for i in range(vsite_count, len(b)):  # This is the last section in the input .itp file
    #     vsite_count += 1

    # while b[vsite_count] != '\n':
    #     vsite_count += 1
    # we need to make a new list of virtual sites that does not include in the sites which have been turned to real atoms
    # term and H_new1 contain the indices of the hydrogens which need to be removed from the virtual sites list

    virtual_sites = []
    if int(args.iteration) == 0:
        for i in range(vsite_index + 2, len(b)):
            site = int(b[i][0:8])
            if site not in H_new2 and site not in H_new1:  # if the site is in the old virtual sites but not in H_new1 or term
                virtual_sites.append(b[i])  # then add it to virtual_sites since those hydrogens are still dummies
    else:
        for i in range(vsite_index + 2, len(b)):
            site = int(b[i][0:8])
            if site not in H_new2 and site not in H_new1 and site not in H_new3:  # don't want virtual sites from H_new3
                virtual_sites.append(b[i])  # then add it to virtual_sites since those hydrogens are still dummies

    # eliminate the old virtual_sites4 list

    del b[vsite_index + 2: len(b)]

    # Insert new virtual sites list

    count = vsite_index + 2
    for i in range(0, len(virtual_sites)):
        b.insert(count, virtual_sites[i])
        count += 1

    # Now write the new topology to a new file

    with open('crosslinked_new.itp', 'w') as target:
        for line in b:
            target.write(line)

    tot_double_bonds = len(c1_atoms)*tot_monomers
    terminated = 0
    for i in range(atoms_index + 2, atoms_count):
        if str.strip(b[i])[-1] == 'T' and str.strip(b[i + 1])[-1] == 'T':
            terminated += 1

    percent_completion = (old_div(float(terminated),float(tot_double_bonds)))*100

    end = time.time()

    # Write important output to a log file

    if int(args.iteration) != 0:
        f = open('xlink.log', 'w')
        f.writelines(['Total time elapsed: %s seconds\n' %(end - start), '\n', 'c1: %s\n' % c1, '\n', 'c2: %s\n' % c2,
                      '\n', 'radical_c2: %s\n' % radical_c2, '\n', 'term: %s\n' % term, '\n', 'term_c1: %s\n' % term_c1,
                      '\n', 'other_c1: %s\n' % other_c1, '\n', 'H_new1: %s\n' % H_new1, '\n', 'H_new2: %s\n' % H_new2,
                      '\n', 'Reacted c2 Radical Indices: %s\n' % c1_rad_ndx, '\n',
                      'c1s which reacted with radical c2: %s\n' % c2_rad_ndx, '\n', 'term_rad_c2: %s\n' % term_rad_c2,
                      '\n', 'reactive_c2_term: %s\n' % reactive_c2_term, '\n', 'H_new3: %s\n' % H_new3, '\n',
                      'Cutoff distance (nm) : %s\n' % cutoff,'\n', 'Vinyl groups terminated: %s\n' % terminated, '\n',
                      'Percent Completion: ' + '{:.1f}'.format(percent_completion) + ' %\n', '\n',
                      'Total Crosslinks: %s\n' % xlinks, 'Stop on Next Iteration?: %s\n' % Stop_next_iter])

    if int(args.iteration) == 0:
        f = open('xlink.log', 'w')
        f.writelines(['Total time elapsed: %s seconds\n' % (end - start), '\n', 'c1: %s\n' % c1, '\n', 'c2: %s\n' % c2,
                      '\n', 'radical_c2: %s\n' % radical_c2, '\n', 'term: %s\n' % term, '\n', 'term_c1: %s\n' % term_c1,
                      '\n', 'other_c1: %s\n' % other_c1, '\n', 'H_new1: %s\n' % H_new1, '\n',
                      'H_new2: %s\n' % H_new2, '\n', 'Cutoff distance (nm) : %s\n' % cutoff, '\n',
                      'Vinyl groups terminated: %s\n' % terminated, '\n',
                      'Percent Completion: ' + '{:.1f}'.format(percent_completion) + '%\n',
                      '\n', 'Total Crosslinks: %s\n' % xlinks, 'Stop on Next Iteration?: %s\n' % Stop_next_iter])
