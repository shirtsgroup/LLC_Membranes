#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import random
from LLC_Membranes.llclib import transform, file_rw, topology
import subprocess
import copy
import os, glob
import sys
import time
import tqdm

location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))  # Directory this script is in


def initialize():

    parser = argparse.ArgumentParser(description='Add water to the tail region of a HII LLC structure')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of structure file where water will be added')
    parser.add_argument('-p', '--top', default='topol.top', help='Name of topology that needs to be updated')
    parser.add_argument('-o', '--out', default='solv_tails.gro', help='Name of output file')
    parser.add_argument('-m', '--monomer', default='NAcarb11V.gro', help='Name of monomer used to build system')
    # ref now added to annotatd .gro file with 'T'
    #parser.add_argument('-r', '--ref', nargs='+', default=['O5', 'O6', 'O7', 'O8', 'O9', 'O10'], help='Names of atoms'
    #                    'where water molecules will be placed in proximity to')
    parser.add_argument('-rmin', default=0.2, type=float, help='Minimum distance away from reference atom to place'
                                                               'water molecules (nm)')
    parser.add_argument('-rmax', default=0.4, type=float, help='Maximum distance away from reference atom to place'
                                                               'water molecules (nm)')
    parser.add_argument('-n', '--nwater', default=600, type=int, help='Number of waters to add to the system')
    parser.add_argument('--restrained', action="store_true", help='Add flag if system contains position restraints')

    args = parser.parse_args()

    return args


def write_em_mdp(steps):
    """
    Write energy minimization .mdp file
    :param steps: number of steps to take using steepest descents algorithm
    :return: Directly writes an energy minimization .mdp file
    """

    with open('em.mdp', 'w') as f:

        f.write("title = Energy Minimization\n")
        f.write("integrator = steep\n")
        f.write("nsteps = %s\n" % steps)
        f.write("cutoff-scheme = verlet\n")
        f.write("nstlist = 40\n")


def random_pt_spherical_shell(pt, rmin, rmax):
    """
    Pick a random point in the region bounded by two spheres
    :param pt: center of concentric spheres
    :param rmin: inside sphere radius
    :param rmax: outside sphere radius
    :return: randomly chosen point between spheres with radii rmin and rmax
    """

    # randomly choose a point within the shell created between two spheres of radius rmin and rmax centered at (0,0,0)
    random_pt = np.random.normal(size=3)  # random vector chosen from gaussian since gaussian is spherically symmetric
    random_pt /= np.linalg.norm(random_pt)  # normalize
    random_pt *= random.uniform(rmin, rmax)  # randomly choose r distance between rmin and rmax

    # translate to random_pt to be centered at pt
    random_pt += pt

    return random_pt

#
# def random_water_orientation(water_xyz, water_alignment_vector, placement):
#     """
#     Randomly orient a water molecule and then place it a desired location
#     :param water_xyz: 3d coordinates of a water molecule
#     :param water_alignment_vector: A reference vector to rotate the water molecule about
#     :param placement: where to place final water configuration in space
#     :return: coordinates of oriented and translated water molecule
#     """
#
#     u = np.random.normal(size=3)  # random vector. From normal distribution since sphere
#     u /= np.linalg.norm(u)  # normalize
#
#     R = transform.Rvect2vect(water_alignment_vector, u)  # rotation matrix to align water_alignment_vector with u
#
#     water_xyz -= water_xyz[0, :]  # center at origin
#
#     rotated = np.zeros([water_xyz.shape[1], 3])
#     for i in range(water_xyz.shape[1]):
#         rotated[i, :] = np.dot(R, water_xyz[i, :])
#
#     rotated += placement  # translate to deisred location
#
#     return rotated


class System(object):

    def __init__(self, gro, top, monomer, rbounds=[0.3, 1], restraints=False, mpi=False, nproc=4):
        """Add set number of water molecules to the tail region of a given configuration

        :param gro: initial coordinate file to be solvated
        :param top: topology associated with gro
        :param nwater: number of water molecules to add to tails
        :param ref: Names of atoms where water molecules will be placed in proximity to
        :param rbounds: list with min, max (in that order) distance of water molecules from a reference atom (nm)

        :type gro: str
        :type top: str
        :type nwater: int
        :type ref: list
        :type rbounds: list

        """

        self.t = md.load(gro)
        self.coordinates = self.t.xyz[0, :, :]  # initial coordinates
        self.natoms = self.t.n_atoms  # number of atoms in the system
        self.rmin, self.rmax = rbounds
        self.restraints = restraints
        self.mpi = mpi
        self.np = nproc

        self.gmx = "gmx"
        if self.mpi:
            self.gmx = "mpirun -np %s gmx_mpi" % self.np

        # handle mdtraj renaming SOL to HOH
        self.res = []
        self.ids = []
        for a in self.t.topology.atoms:
            if a.residue.name == 'HOH':
                self.res.append('SOL')
                if a.name == 'H1':
                    self.ids.append('HW1')
                elif a.name == 'H2':
                    self.ids.append('HW2')
                elif a.name == 'O':
                    self.ids.append('OW')
            else:
                self.res.append(a.residue.name)
                self.ids.append(a.name)

        self.full_box = self.t.unitcell_vectors  # unitcell vectors
        self.box_gromacs = [self.full_box[0, 0, 0], self.full_box[0, 1, 1], self.full_box[0, 2, 2],
                            self.full_box[0, 0, 1], self.full_box[0, 2, 0], self.full_box[0, 1, 0],
                            self.full_box[0, 0, 2], self.full_box[0, 1, 2], self.full_box[0, 2, 0]]

        self.water = topology.Molecule('HOH').xyz[0, ...]
        self.water_alignment_vector = self.water[0, :] - np.mean(self.water, axis=0)  # vector around which water molecule can be rotated

        self.topname = top

        with open(self.topname, 'r') as f:
            self.top = []
            for line in f:
                self.top.append(line)

        self.add_water_placeholder()

        self.monomer = topology.LC(monomer)
        self.ref_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.tail_atoms]  # indices of reference atoms
        self.ref_atom_locations = self.t.xyz[0, self.ref_atoms, :]  # coordinates of reference atoms

        self.placement_options = [i for i in range(len(self.ref_atoms))]  # list of indices corresponding to atoms in ref_atom_locations

        # need to remove hardcoding. Probably use llclib.topology.Residue
        self.water_ids = ['OW', 'HW1', 'HW2']  # water atom names
        self.water_res = ['SOL', 'SOL', 'SOL']  # water residue name

        write_em_mdp(5)  # 5 step energy minimization

    def add_water_placeholder(self):
        """
        Open up topology, add water include statement and a PLACEHOLDER so that varying amounts of water can be placed
        :param top: Name of topology to be modified
        :return: placeholder_top.top is written
        """

        solvated = False
        for i in range(len(self.top)):
            if self.top[i].count('Water Topology') >= 1 or self.top[i].count('SOL Topology') >= 1:
                solvated = True

        if not solvated:
            # find [ system ] directive
            system_ndx = 0
            while self.top[system_ndx].count('[ system ]') == 0:
                system_ndx += 1

            self.top.insert(system_ndx, '; Water Topology\n')
            self.top.insert(system_ndx + 1, '#include "%s/../top/Forcefields/gaff/tip3p.itp"\n' % location)
            self.top.insert(system_ndx + 2, '\n')

        self.top.append('SOL                PLACEHOLDER')

        with open('placeholder_top.top', 'w') as f:

            for line in self.top:
                f.write(line)

    def place_water(self, nwater, steps):
        """
        Place water molecules near a random reference atom and perform a short energy minimization
        :param placement_options: All of the possible indexes of reference atoms near which a water molecule could be placed (list)
        :param ref_atom_locations: xyz coordinates of reference atoms (numpy array [number of ref atoms, 3]
        :param coordinates: Coordinates of full system (numpy array [natoms, 3])
        :param rmin: minimimum distance to place water molecules from reference atom (float)
        :param rmax: maximum distance to place water molecules from reference atom (float)
        :param ids: all atom names in order they appear in coordinates (list)
        :param res: all residue atoms in order they appear in coordinates (list)
        :param nwater: number of water molecules in the system (int)
        :param steps: number of energy minimization steps to take (int)
        :return: Potential energy of slightly minimized system, coordinates of minimized system, atom used for placement,
        new locations of reference atoms
        """

        placement_atom = np.random.choice(self.placement_options)  # choose which atom to place water molecule near
        placement = random_pt_spherical_shell(self.ref_atom_locations[placement_atom, :], self.rmin,
                                                          self.rmax)  # point near placement atom
        placed_water_coordinates = transform.random_orientation(self.water, self.water_alignment_vector, placement)
        new_coordinates = np.concatenate((self.coordinates, placed_water_coordinates),
                                         axis=0)  # add to full list of coordinates
        names = self.ids + self.water_ids  # add water to ids
        residues = self.res + self.water_res  # add water residue to res
        file_rw.write_gro_pos(new_coordinates, 'water.gro', ids=names, res=residues,
                              box=self.box_gromacs)  # write out config with new water
        minimized_coordinates, new_ref_atom_locations = self.energy_minimize(steps, nwater + 1)  # energy minimzed system
        nrg = subprocess.check_output(
            ["awk", "/Potential Energy/ {print $4}", "em.log"])  # get Potential energy from em.log

        return float(nrg.decode("utf-8")), minimized_coordinates, placement_atom, new_ref_atom_locations

    def energy_minimize(self, steps, nwater):
        """
        Energy minimize a configuration
        :param steps: number of steepest descent energy minimization steps to take
        :param nwater: number of water molecules in the system
        :return: coordinates of energy minimized structure, updated coordinates of reference atoms
        """

        write_em_mdp(steps)  # write em.mdp with a given number of steps

        cp = "cp placeholder_top.top top_intermediate.top"
        p1 = subprocess.Popen(cp.split())  # make a copy of placeholder_top.top
        p1.wait()

        sed = "sed -i -e s/PLACEHOLDER/%s/g top_intermediate.top" % nwater
        p2 = subprocess.Popen(sed.split())
        p2.wait()

        grompp = "%s grompp -f em.mdp -p top_intermediate.top -c water.gro -o em" % self.gmx
        if self.restraints:
            grompp += " -r water.gro"
        p3 = subprocess.Popen(grompp.split(), stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        p3.wait()

        mdrun = "%s mdrun -deffnm em" % self.gmx
        p4 = subprocess.Popen(mdrun.split(), stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        p4.wait()

        t = md.load('em.gro')
        minimized_coordinates = t.xyz[0, :, :]  # coordinates of energy minimized system
        new_ref_atom_locations = t.xyz[0, self.ref_atoms, :]  # new coordinates of reference atoms

        return minimized_coordinates, new_ref_atom_locations

    def insert_all_water(self, nwater, output='solvated.gro', final_topname='topol.top'):

        trials = 0  # number of water placement tries
        # start = time.time()

        for i in tqdm.tqdm(range(nwater), unit='waters'):
            # randomly place water molecule and do short energy minimization
            energy, new_coordinates, placement_atom, new_ref_atom_locations = self.place_water(i, 5)
            trials += 1
            count = 0
            while energy > -50000:  # make sure the potential energy doesn't get too close to exploding

                energy, new_coordinates, placement_atom, new_ref_atom_locations = self.place_water(i, 5)
                trials += 1
                count += 1
                if count > 10:
                    # if the energy is too high for water placement, run a longer energy minimization
                    sys.stdout.write("\r Running longer energy minimization...                                         "
                                     "                \r")
                    sys.stdout.flush()
                    file_rw.write_gro_pos(self.coordinates, 'water.gro', box=self.box_gromacs, ids=self.ids, res=self.res)
                    self.coordinates, self.ref_atom_locations = self.energy_minimize(500, i)  # energy minimize and update coordinates
                    count = 0
                    write_em_mdp(1)

            # success_rate = 100*((i + 1) / trials)
            # s = "Waters placed: %s/%s, Placement Success Rate : %3.2f, Potential Energy = %s\r" % (i, nwater,
            #                                                                                        success_rate, energy)
            # sys.stdout.write("\r"+s)
            # sys.stdout.flush()
            self.ids += self.water_ids
            self.res += self.water_res
            self.coordinates = copy.deepcopy(new_coordinates)
            self.placement_options.remove(placement_atom)
            for filename in glob.glob("./#*"):
                os.remove(filename)

        cp = "cp top_intermediate.top %s" % final_topname
        p = subprocess.Popen(cp.split())
        p.wait()

        # energy minimize final configuration until convergence
        write_em_mdp(-1)
        grompp = "%s grompp -f em.mdp -p %s -c water.gro -o %s" % (self.gmx, final_topname, output.split('.')[0])
        if self.restraints:
            grompp += " -r water.gro"
        p = subprocess.Popen(grompp.split())
        p.wait()

        mdrun = "%s mdrun -v -deffnm %s" % (self.gmx, output.split('.')[0])
        p = subprocess.Popen(mdrun.split())
        p.wait()


if __name__ == "__main__":

    args = initialize()

    if os.path.exists(args.out):  # Get rid of output file if it already exists since it can mess things up
        os.remove(args.out)

    system = System(args.gro, args.top, args.monomer, rbounds=[args.rmin, args.rmax], restraints=args.restrained)
    system.insert_all_water(args.nwater)

    # water = md.load('%s/../top/solutes/water.gro')  # load water structure
    # water_xyz = water.xyz[0, :, :]  # get water coordinates
    # water_centroid = np.mean(water_xyz, axis=0)  # get centroid of water
    # water_alignment_vector = water_xyz[0, :] - water_centroid  # vector around which water molecule can be rotated

    # add_water_placeholder(top)  # add a placeholder field to topol.top

    # ref_atoms = [a.index for a in t.topology.atoms if a.name in args.ref]  # indices of reference atoms
    #
    # ref_atom_locations = t.xyz[0, ref_atoms, :]  # coordinates of reference atoms
    #
    # placement_options = [i for i in range(len(ref_atoms))]  # list of indices corresponding to atoms in ref_atom_locations

    # rmin = args.rmin  # min distance from ref atom to place water molecule
    # rmax = args.rmax  # max distance from ref atom to place water molecule


    # trials = 0  # number of water placement tries
    #
    # start = time.time()  # when water placement starts
    #
    # for i in range(nwater):
    #     # randomly place water molecule and do short energy minimization
    #     energy, new_coordinates, placement_atom, new_ref_atom_locations = place_water(placement_options,
    #                                                     ref_atom_locations, coordinates, rmin,  rmax, ids, res, i, 5)
    #     trials += 1
    #     count = 0
    #     while energy > -50000:  # make sure the potential energy doesn't get too close to exploding
    #
    #         energy, new_coordinates, placement_atom, new_ref_atom_locations = place_water(placement_options,
    #                                                         ref_atom_locations, coordinates, rmin, rmax, ids, res, i, 5)
    #         trials += 1
    #         count += 1
    #         if count > 10:
    #             # if the energy is too high for water placement, run a longer energy minimization
    #             sys.stdout.write("\r Running longer energy minimization...                                             "
    #                              "             \r")
    #             sys.stdout.flush()
    #             file_rw.write_gro_pos(coordinates, 'water.gro', box=box_gromacs, ids=ids, res=res)
    #             coordinates, ref_atom_locations = energy_minimize(500, i)  # energy minimize and update coordinates
    #             count = 0
    #             write_em_mdp(1)
    #
    #     success_rate = 100*((i + 1) / trials)
    #     s = "Waters placed: %s/%s, Placement Success Rate : %3.2f, Potential Energy = %s\r" % (i, nwater, success_rate, energy)
    #     sys.stdout.write("\r"+s)
    #     sys.stdout.flush()
    #     ids += water_ids
    #     res += water_res
    #     coordinates = copy.deepcopy(new_coordinates)
    #     placement_options.remove(placement_atom)
    #     for filename in glob.glob("./#*"):
    #         os.remove(filename)

    # print("Waters placed in %4.2f seconds" % (time.time() - start))
    # file_rw.write_gro_pos(coordinates, args.out, box=box_gromacs, res=res, ids=ids)
