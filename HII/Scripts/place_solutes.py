#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
import Atom_props
from llclib import file_rw
from llclib import transform
import subprocess
import os
from gentop import SystemTopology
import tqdm


script_location = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


def initialize():

    parser = argparse.ArgumentParser(description='Add specified amount of solvent to box')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file to add solutes to')
    parser.add_argument('-c', '--concentration', nargs='+', type=float, help='Concentration of solute (M)')
    parser.add_argument('-n', '--n_solute', nargs='+', help='Number of solute molecules to add (overrides '
                                                            'concentration')
    parser.add_argument('-s', '--solutes', nargs='+', help='.gro file for solute molecules')

    args = parser.parse_args()

    return args


def concentration_to_nsolute(conc, box_vectors, solute):
    """
    :param conc: (float) desired solute concentration (M)
    :param box_vectors: (numpy array, (3, 3)) box vectors. Each row represents a box vector.
    :param solute: mdtraj trajectory object generated from solute configuration file (.gro)
    :return: (int) number of solute molecules to add to box to achieve desired concentration
    """

    V = np.dot(box_vectors[2, :], np.cross(box_vectors[0, :], box_vectors[1, :]))  # box volume (nm^3)
    V *= 1 * 10 ** -24  # convert to L
    mols_solute = conc * V  # number of mols of solvent to add

    # mw = 0  # molecular weight (grams)
    # for a in solute.topology.atoms:
    #     mw += Atom_props.mass[a.name]

    mass_to_add = solute.mw * mols_solute

    NA = 6.022 * 10 ** 23  # avogadro's number
    mass_solute = solute.mw / NA  # mass of a single solutes (grams)

    nsolute = int(mass_to_add / mass_solute)  # number of solute molecules to add

    actual_concentration = nsolute / (NA*V)  # mol/L

    return nsolute, actual_concentration


def random_orientation(xyz, alignment_vector, placement):
    """
    Randomly orient a water molecule and then place it a desired location
    :param water_xyz: 3d coordinates of a water molecule
    :param water_alignment_vector: A reference vector to rotate the water molecule about
    :param placement: where to place final water configuration in space
    :return: coordinates of oriented and translated water molecule
    """

    u = np.random.normal(size=3)  # random vector. From normal distribution since sphere
    u /= np.linalg.norm(u)  # normalize

    R = transform.Rvect2vect(alignment_vector, u)  # rotation matrix to align water_alignment_vector with u

    xyz -= xyz[0, :]  # center at origin

    rotated = np.zeros([xyz.shape[1], 3])
    for i in range(xyz.shape[1]):
        rotated[i, :] = np.dot(R, xyz[i, :])

    rotated += placement  # translate to deisred location

    return rotated


def net_charge(nsolute, solutes):
    """
    :param nsolute: list of number of solutes to be added
    :param solutes: list of solute objects
    :return: net charge of system after addition of nsolute
    """

    net_charge = 0
    for i, n in enumerate(nsolute):
        net_charge += n*solutes[i].charge

    return net_charge


class Solvent(object):

    def __init__(self, gro, intermediate_fname='solvate.gro', em_steps=100, p_coupling='isotropic'):
        """
        :param gro: configuration of solvent
        :param intermediate_fname : name of intermediate .gro files if placing solute in box
        :param em_steps : number of energy minimization steps if placing solute in box
        """

        t = md.load(gro)
        self.box_vectors = t.unitcell_vectors[0, :, :]  # box vectors

        self.box_gromacs = [self.box_vectors[0, 0], self.box_vectors[1, 1], self.box_vectors[2, 2],
                            self.box_vectors[0, 1], self.box_vectors[2, 0], self.box_vectors[1, 0],
                            self.box_vectors[0, 2], self.box_vectors[1, 2], self.box_vectors[2, 0]]  # box in gromacs format

        self.positions = t.xyz[0, :, :]  # positions of all atoms
        self.residues = []
        self.names = []
        self.top = SystemTopology(gro)
        self.intermediate_fname = intermediate_fname
        self.em_steps = em_steps

        # because mdtraj changes the names
        for a in t.topology.atoms:
            if a.residue.name == 'HOH':
                self.residues.append('SOL')
                if a.name == 'O':
                    self.names.append('OW')
                elif a.name == 'H1':
                    self.names.append('HW1')
                elif a.name == 'H2':
                    self.names.append('HW2')
            else:
                self.residues.append(a.residue.name)
                self.names.append(a.name)

    def place_solute(self, solute):
        """
        :param solute: Solute object generated from solute configuration file (.gro)
        """
        placement_point = self.random_point_box()  # where to place solute
        solute_positions = transform.translate(solute.xyz[0, :, :], solute.xyz[0, 0, :], placement_point)  # translate solute to placement point
        self.positions = np.concatenate((self.positions, solute_positions))  # add to array of positions
        self.residues += solute.residues  # add solute residues to list of all residues
        self.names += solute.names  # add solute atom names to list of all names
        self.top.add_residue(solute, write=True)  # add 1 solute to topology

        # write new .gro file
        file_rw.write_gro_pos(self.positions, self.intermediate_fname, box=self.box_gromacs, ids=self.names, res=self.residues)
        nrg = self.energy_minimize(self.em_steps)

        if nrg >= 0:
            self.revert(solute)
            self.place_solute(solute)

    def energy_minimize(self, steps):
        """
        Energy minimize a configuration
        :param steps: number of steepest descent energy minimization steps to take
        :return: coordinates of energy minimized structure, updated coordinates of reference atoms
        """

        file_rw.write_em_mdp(steps)  # write em.mdp with a given number of steps

        p1 = subprocess.Popen(
            ["gmx", "grompp", "-p", "topol.top", "-f", "em.mdp", "-o", "em", "-c", "%s" % self.intermediate_fname],
            stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # generate atomic level input file
        p1.wait()
        p2 = subprocess.Popen(["gmx", "mdrun", "-deffnm", "em"], stdout=open(os.devnull, 'w'),
                              stderr=subprocess.STDOUT)  # run energy minimization
        p2.wait()

        nrg = subprocess.check_output(
            ["awk", "/Potential Energy/ {print $4}", "em.log"])  # get Potential energy from em.log

        try:
            return float(nrg.decode("utf-8"))
        except ValueError:
            return 0  # If the system did not energy minimize, the above statement will not work because nrg will be an
            # empty string. Make nrg=0 so placement gets attempted again

    def random_point_box(self):
        """
        :param box_vectors: (numpy array, (3, 3)) box vectors. Each row represents a box vector.
        :return: (numpy array, (3)) coordinates of a randomly chosen point that lies in box
        """

        A = self.box_vectors[0, :]  # x box vector
        B = self.box_vectors[1, :]  # y box vector
        C = self.box_vectors[2, :]  # z box vector
        u, v, w = np.random.rand(3)  # generate 3 random numbers between 0 and 1
        pt = np.array([0, 0, 0]) + u * A + v * B + w * C  # places point inside 3D box defined by box vector A, B and C

        return pt

    def revert(self, solute):
        """
        Revert system to how it was before solute addition
        :return:
        """
        n = -solute.natoms
        self.positions = self.positions[:n, :]
        self.residues = self.residues[:n]
        self.names = self.names[:n]
        self.top.add_residue(solute, n=-1, write=False)  # subtract a solute from the topology


class Solute(object):

    def __init__(self, name):

        self.is_ion = False
        # check if residue is an ion
        with open('%s/../top/topologies/ions.txt' % script_location) as f:
            ions = []
            for line in f:
                if line[0] != '#':
                    ions.append(str.strip(line))

        if name in ions:
            self.is_ion = True
            self.residues = [name]
            self.names = [name]
            self.xyz = np.zeros([1, 1, 3])
            self.xyz[0, 0, :] = [0, 0, 0]
            self.natoms = 1
            self.mw = Atom_props.mass[name]
            self.charge = Atom_props.charge[name]
            self.resname = name
        else:
            try:
                t = md.load('%s.pdb' % name, standard_names=False)  # see if there is a solute configuration in this directory
            except OSError:
                try:
                    t = md.load('%s/../top/topologies/%s.pdb' % (script_location, name), standard_names=False)  # check if the configuration is
                    # located with all of the other topologies
                except OSError:
                    print('No residue %s found' % name)
                    exit()

            try:
                f = open('%s.itp' % name, 'r')
            except FileNotFoundError:
                try:
                    f = open('%s/../top/topologies/%s.itp' % (script_location, name), 'r')
                except FileNotFoundError:
                    print('No topology %s.itp found' % name)

            itp = []
            for line in f:
                itp.append(line)

            f.close()

            self.natoms = t.n_atoms

            atoms_index = 0
            while itp[atoms_index].count('[ atoms ]') == 0:
                atoms_index += 1

            atoms_index += 2
            self.charge = 0
            for i in range(self.natoms):
                self.charge += float(itp[atoms_index + i].split()[6])

            self.residues = [a.residue.name for a in t.topology.atoms]
            self.resname = self.residues[0]
            self.names = [a.name for a in t.topology.atoms]
            self.xyz = t.xyz

            self.mw = 0  # molecular weight (grams)
            for a in t.topology.atoms:
                self.mw += Atom_props.mass[a.name]


if __name__ == "__main__":

    args = initialize()

    os.environ["GMX_MAXBACKUP"] = "-1"  # stop GROMACS from making backups

    solvent = Solvent(args.gro)
    if args.concentration:
        concentration = [float(i) for i in args.concentration]
    elif args.n_solute:
        n = [int(i) for i in args.n_solute]

    solutes = []
    nsolute = []

    for i, s in enumerate(args.solutes):

        solutes.append(Solute(s))

        if args.concentration:
            n, actual_concentration = concentration_to_nsolute(concentration[i], solvent.box_vectors, solutes[i])
            nsolute.append(n)
            print("Actual Concentration of %s : %.2f mol/L" % (s, actual_concentration))
        elif args.n_solute:
            nsolute[i] = n[i]
        else:
            print("You must specify a concentration or number of solute molecules")
            exit()

    # system_charge = net_charge(nsolute, solutes)
    # solute_charges = np.array([int(a.charge) for a in solutes])
    #
    # # ensure charge neutrality
    # while -0.0001 > system_charge or system_charge > 0.0001:
    #     if system_charge < 0:
    #         nsolute[np.random.choice(np.where(solute_charges > 0.0)[0])] += 1
    #     elif system_charge > 0:
    #         nsolute[np.random.choice(np.where(solute_charges < 0.0)[0])] += 1
    #     system_charge = net_charge(nsolute, solutes)

    # a "smarter" way to add solutes that keeps things relatively neutral - has bugs
    # while sum(nsolute) > 0:
    #     print(sum(nsolute))
    #     solute_to_add = np.random.choice(np.where(nsolute != 0)[0])  # randomly choose a solute to add
    #     nsolute[solute_to_add] -= 1  # subtract from solutes that need to be added
    #     solvent.place_solute(solutes[solute_to_add])  # place the solute
    #     if solutes[solute_to_add].charge > 0:  # if it is a cation, neutralize with anion(s)
    #         # choose anion such that its charge is positive but not greater in magnitude than the cation already added
    #         anion = np.random.choice(np.where(solute_charges < 0.0)[0])
    #         nadd = int(solutes[anion].charge / solutes[solute_to_add].charge)  # if charge on cation is +2 and charge on anion is -1, add two anions
    #         for i in range(nadd):
    #             solvent.place_solute(solutes[anion])
    #             nsolute[anion] -= 1
    #     elif solutes[solute_to_add].charge < 0:  # if it is an anion, neutralize with cation(s)
    #         # choose cation such that its charge is positive but not greater in magnitude than the anion already added
    #         cation = np.random.choice(np.where(solute_charges > 0.0)[0])
    #         nadd = int(solutes[cation].charge / solutes[solute_to_add].charge)  # if charge on cation is +2 and charge on anion is -1, add two anions
    #         for i in range(nadd):
    #             solvent.place_solute(solutes[cation])
    #             nsolute[cation] -= 1

    # print(int(solutes[i].charge))
    # print("Adding %d %s molecules" % (nsolute, s))

    for n in tqdm.tqdm(range(len(nsolute))):
        solvent.place_solute(solutes[n])

    from pathlib import Path

    for p in Path(".").glob("step*"):
        p.unlink()