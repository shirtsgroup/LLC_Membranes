#!/usr/bin/env python

import argparse
import mdtraj as md
import numpy as np
from LLC_Membranes.analysis import Atom_props
from LLC_Membranes.llclib import file_rw, transform
from LLC_Membranes.setup.gentop import SystemTopology
import subprocess
import os
import tqdm
import matplotlib.path as path
from scipy import spatial

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

    rotated += placement  # translate to desired location

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


def put_in_box(pt, x_box, y_box, m, angle):
    """
    :param pt: The point to place back in the box
    :param x_box: length of box in x dimension
    :param y_box: length of box in y dimension
    :param m: slope of box vector
    :param angle: angle between x axis and y box vector
    :return: coordinate shifted into box
    """

    b = - m * x_box  # y intercept of box vector that does not pass through origin (right side of box)
    if pt[1] < 0:
        pt[:2] += [np.cos(angle)*x_box, np.sin(angle)*x_box]  # if the point is under the box
    if pt[1] > y_box:
        pt[:2] -= [np.cos(angle)*x_box, np.sin(angle)*x_box]
    if pt[1] > m*pt[0]:  # if the point is on the left side of the box
        pt[0] += x_box
    if pt[1] < m*(pt[0] - b):  # if the point is on the right side of the box
        pt[0] -= x_box

    return pt


def trace_pores(pos, box, layers):
    """
    Find the line which traces through the center of the pores
    :param pos: positions of atoms used to define pore location (args.ref) [natoms, 3]
    :param box: xy box vectors, [2, 2], mdtraj format
    :param layers: number of layers
    :return: points which trace the pore center
    """

    atoms_p_pore = int(pos.shape[0] / 4)  # atoms in each pore
    atoms_p_layer = int(atoms_p_pore / layers)  # atom per layer

    v = np.zeros([4, 2])  # vertices of unitcell box
    v[0, :] = [0, 0]
    v[1, :] = [box[0, 0], 0]
    v[3, :] = [box[1, 0], box[1, 1]]
    v[2, :] = v[3, :] + [box[0, 0], 0]

    center = [np.mean(v[:, 0]), np.mean(v[:, 1]), 0]  # geometric center of box
    bounds = path.Path(v)  # create a path tracing the vertices, v

    angle = np.arccos(box[1, 1]/box[0, 0])  # angle of monoclinic box
    if box[1, 0] < 0:  # the case of an obtuse angle
        angle += np.pi / 2

    m = (v[3, 1] - v[0, 1]) / (v[3, 0] - v[0, 0])  # slope from points connecting first and third vertices

    centers = np.zeros([4*layers, 3])

    for p in range(4):
        pore = pos[p*atoms_p_pore:(p+1)*atoms_p_pore, :]  # coordinates for atoms belonging to a single pore
        for l in range(layers):
            before = pore[l*atoms_p_layer, :]  # choose the first atom as a reference

            shift = transform.translate(pore[l*atoms_p_layer:(l+1)*atoms_p_layer, :], before, center)  # shift everything to towards the center

            for i in range(shift.shape[0]):  # check if the points are within the bounds of the unitcell
                if not bounds.contains_point(shift[i, :2]):
                    shift[i, :] = put_in_box(shift[i, :], box[0, 0], box[1, 1], m, angle)  # if its not in the unitcell, shift it so it is

            c = np.zeros([1, 3])
            c[0, :] = [np.mean(shift[:, 0]), np.mean(shift[:, 1]), np.mean(shift[:, 2])]  # geometric center of reference atoms in this layer

            centers[p*layers + l, :] = transform.translate(c, center, before)  # move everything back to where it was

            if not bounds.contains_point(centers[p*layers, :]):  # make sure everything is in the box again
                centers[p*layers + l, :] = put_in_box(centers[p*layers + l, :], box[0, 0], box[1, 1], m, angle)

    return centers


def placement(z, pts, box):
    """
    :param z: z location where solute should be placed
    :param pts: points which run through the pore
    :return: location to place solute
    """

    # check if point is already in the spline
    if z in pts[:, 2]:

        ndx = np.where(pts[:, 2] == z)[0][0]

        return pts[ndx, :]
    # otherwise interpolate between closest spline points
    else:
        v = np.zeros([4, 2])  # vertices of unitcell box
        v[0, :] = [0, 0]
        v[1, :] = [box[0, 0], 0]
        v[3, :] = [box[1, 0], box[1, 1]]
        v[2, :] = v[3, :] + [box[0, 0], 0]
        center = [np.mean(v[:, 0]), np.mean(v[:, 1]), 0]  # geometric center of box

        bounds = path.Path(v)  # create a path tracing the vertices, v

        angle = np.arccos(box[1, 1]/box[0, 0])  # angle of monoclinic box
        if box[1, 0] < 0:  # the case of an obtuse angle
            angle += np.pi / 2

        m = (v[3, 1] - v[0, 1]) / (v[3, 0] - v[0, 0])  # slope from points connecting first and fourth vertices

        # shift = transform.translate(z, before, center)
        #
        # put_in_box(pt, box[0, 0], box[1, 1], m, angle)

        # find z positions, in between which solute will be placed
        lower = 0
        while pts[lower, 2] < z:
            lower += 1

        upper = pts.shape[0] - 1
        while pts[upper, 2] > z:
            upper -= 1

        limits = np.zeros([2, 3])
        limits[0, :] = pts[lower, :]
        limits[1, :] = pts[upper, :]

        shift = transform.translate(limits, limits[0, :], center)  # shift limits to geometric center of unit cell
        shift[:, 2] = [limits[0, 2], limits[1, 2]]  # keep z positions the same

        for i in range(shift.shape[0]):  # check if the points are within the bounds of the unitcell
            if not bounds.contains_point(shift[i, :2]):
                shift[i, :] = put_in_box(shift[i, :], box[0, 0], box[1, 1], m, angle)

        # Use parametric representation of line between upper and lower points to find the xy value where z is satsified
        v = shift[1, :] - shift[0, :]  # direction vector

        t = (z - shift[0, 2]) / v[2]  # solve for t since we know z
        x = shift[0, 0] + t*v[0]
        y = shift[0, 1] + t*v[1]

        place = np.zeros([1, 3])
        place[0, :] = [x, y, 0]
        place = transform.translate(place, center, limits[0, :])  # put xy coordinate back
        place[0, 2] = z

        if not bounds.contains_point(place[0, :]):  # make sure everything is in the box again
            place[0, :] = put_in_box(place[0, :], box[0, 0], box[1, 1], m, angle)

        return place[0, :]


class Solvent(object):

    def __init__(self, gro, intermediate_fname='solvate.gro', em_steps=100, p_coupling='isotropic'):
        """
        :param gro: configuration of solvent
        :param intermediate_fname : name of intermediate .gro files if placing solute in box
        :param em_steps : number of energy minimization steps if placing solute in box
        """

        self.t = md.load(gro)
        self.box_vectors = self.t.unitcell_vectors[0, :, :]  # box vectors

        self.box_gromacs = [self.box_vectors[0, 0], self.box_vectors[1, 1], self.box_vectors[2, 2],
                            self.box_vectors[0, 1], self.box_vectors[2, 0], self.box_vectors[1, 0],
                            self.box_vectors[0, 2], self.box_vectors[1, 2], self.box_vectors[2, 0]]  # box in gromacs format

        self.positions = self.t.xyz[0, :, :]  # positions of all atoms
        self.residues = []
        self.names = []
        self.top = SystemTopology(gro)
        self.intermediate_fname = intermediate_fname
        self.em_steps = em_steps

        # data specifically required for adding solutes to pores
        self.pore_spline = None
        self.water = [a.index for a in self.t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']
        self.water_top = Solute('SOL')

        # because mdtraj changes the names
        for a in self.t.topology.atoms:
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

    def place_solute(self, solute, placement_point, random=False, freeze=False, rem=.5):

        """
        Place solute at desired point and energy minimze the system
        :param solute: name of solute object (str)
        :param placement_point: point to place solute (np.array([3])
        :param random: place solute at random point in box (bool)
        :param freeze: freeze all atoms outside rem during energy minimization (bool)
        :param rem: radius from placement_point within which atoms will NOT be frozen (float, nm)
        :return:
        """

        # p = placement_point - (solute.com - solute.xyz[0, 0, :])  # want to move com to placement point
        #solute_positions = random_orientation(solute.xyz[0, ...], solute.xyz[0, 0, :] - solute.xyz[0, 1, :], placement_point)  # to be fixed in THE FUTURE
        solute_positions = transform.translate(solute.xyz[0, :, :], solute.com, placement_point)  # translate solute to placement point
        self.positions = np.concatenate((self.positions, solute_positions))  # add to array of positions
        self.residues += solute.residues  # add solute residues to list of all residues
        self.names += solute.names  # add solute atom names to list of all names
        self.top.add_residue(solute, write=True)  # add 1 solute to topology

        # write new .gro file
        file_rw.write_gro_pos(self.positions, self.intermediate_fname, box=self.box_gromacs, ids=self.names, res=self.residues)

        if freeze:
            self.freeze_ndx(solute_placement_point=placement_point, res=solute.resname)

        nrg = self.energy_minimize(self.em_steps, freeze=freeze)

        if nrg >= 0:
            self.revert(solute)
            if random:
                self.place_solute_random(solute)
            else:
                #self.remove_water(placement_point, 3)
                self.place_solute(solute, placement_point, freeze=False)

    def place_solute_random(self, solute):
        """
        :param solute: Solute object generated from solute configuration file (.gro)
        """
        placement_point = self.random_point_box()  # where to place solute
        self.place_solute(solute, placement_point, random=True)

    def place_solute_pores(self, solute, z=None, layers=20, pores=4, ref=['C', 'C1', 'C2', 'C3', 'C4', 'C5']):
        """
        Place solute in middle of pores at given z location
        :param solute: solute object
        :param z: z location of solute center of mass (float)
        :param layers: number of layers in system (when initial configuration was set up) (int)
        :param pores: number of pores in which to place solutes (int)
        :param ref: reference atoms used to define pore center
        :return:
        """

        if not self.pore_spline:
            ref = [a.index for a in self.t.topology.atoms if a.name in ref]
        # redo each time because positions change slightly upon energy minimization
        self.pore_spline = trace_pores(self.positions[ref, :], self.box_vectors[:2, :2], layers)

        # format z so that it is an array
        if type(z) is float or type(z) is np.float64:
            z = np.array([z for i in range(pores)])

        for i in tqdm.tqdm(range(pores)):
            placement_point = placement(z[i], self.pore_spline[i * layers: (i + 1) * layers, :], self.box_vectors[:2, :2])
            self.place_solute(solute, placement_point, freeze=True)

    def energy_minimize(self, steps, freeze=False, freeze_group='Freeze', freeze_dim='xyz'):
        """
        Energy minimize a configuration
        :param steps: number of steepest descent energy minimization steps to take
        :return: coordinates of energy minimized structure, updated coordinates of reference atoms
        """

        # write em.mdp with a given number of steps
        file_rw.write_em_mdp(steps, freeze=freeze, freeze_group='Freeze', freeze_dim='xyz')

        if freeze:
            p1 = subprocess.Popen(
                ["gmx", "grompp", "-p", "topol.top", "-f", "em.mdp", "-o", "em", "-c", "%s" % self.intermediate_fname,
                 "-n", "freeze_index.ndx"],
                stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # generate atomic level input file
            p1.wait()
        else:
            p1 = subprocess.Popen(
                ["gmx", "grompp", "-p", "topol.top", "-f", "em.mdp", "-o", "em", "-c", "%s" % self.intermediate_fname],
                stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)  # generate atomic level input file
            p1.wait()

        p2 = subprocess.Popen(["gmx", "mdrun", "-deffnm", "em"], stdout=open(os.devnull, 'w'),
                              stderr=subprocess.STDOUT)  # run energy minimization
        p2.wait()

        p3 = subprocess.Popen(["cp", "em.gro", "%s" % self.intermediate_fname])
        p3.wait()

        # update positions
        self.positions = md.load('%s' % self.intermediate_fname).xyz[0, :, :]

        nrg = subprocess.check_output(
            ["awk", "/Potential Energy/ {print $4}", "em.log"])  # get Potential energy from em.log

        try:
            return float(nrg.decode("utf-8"))
        except ValueError:
            return 0  # If the system did not energy minimize, the above statement will not work because nrg will be an
            # empty string. Make nrg=0 so placement gets attempted again

    def freeze_ndx(self, solute_placement_point=None, rem=None, res=None):
        """
        Write an index file for atoms to be frozen
        :param solute_placement_point: xyz position of where water molecule was placed
        :param rem: spherical radius measured from water molecule placement point outside which all atoms will be frozen
        :param res: freeze this residue and no other atoms (can be combined with rem option)
        :return: index file with indices of atoms to be frozen
        """

        freeze_indices = []
        if rem:
            pts = spatial.cKDTree(self.positions).query_ball_point(solute_placement_point, rem)
            freeze_indices = [a.index for a in self.t.topology.atoms if a.index not in pts]
        elif res:
            freeze_indices += [a for a in range(len(self.residues)) if self.residues[a] == res]
        else:
            print('WARNING: No valid options supplied in order to determine freeze indices. Specify rem or res.')

        with open('freeze_index.ndx', 'w') as f:

            f.write('[ Freeze ]\n')
            for i, entry in enumerate(freeze_indices):
                if (i + 1) % 15 == 0:
                    f.write('{:5d}\n'.format(entry + 1))
                else:
                    f.write('{:5d} '.format(entry + 1))

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
        """
        n = -solute.natoms
        self.positions = self.positions[:n, :]
        self.residues = self.residues[:n]
        self.names = self.names[:n]
        self.top.add_residue(solute, n=-1, write=False)  # subtract a solute from the topology

    def write_config(self, name='out.gro'):
        """
        Write .gro coordinate file from current positions
        :param name: name of coordinate file to write (str)
        """
        # write new .gro file
        file_rw.write_gro_pos(self.positions, name, box=self.box_gromacs, ids=self.names, res=self.residues)

    def remove_water(self, point, n):

        """
        remove n water molecules closest to point
        """

        tree = spatial.cKDTree(self.positions[self.water, :])
        rm = []

        nn = tree.query(point, k=n)[1]
        for j in nn:
            rm.append(self.water[j])
            rm.append(self.water[j] + 1)
            rm.append(self.water[j] + 2)

        # update relevant arrays
        self.positions = np.delete(self.positions, rm, axis=0)
        self.residues = [self.residues[x] for x in range(len(self.residues)) if x not in rm]
        self.names = [self.names[x] for x in range(len(self.names)) if x not in rm]
        self.water = [i for i, x in enumerate(self.residues) if x == 'SOL' and self.names[i] == 'OW']

        self.top.remove_residue(self.water_top, n, write=True)


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

            self.com = np.zeros([3])  # center of mass of solute
            for i in range(self.xyz.shape[1]):
                self.com += self.xyz[0, i, :] * Atom_props.mass[self.names[i]]
            self.com /= self.mw


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
            nsolute.append(n[i])
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

    # for n in tqdm.tqdm(range(len(nsolute))):
    #     solvent.place_solute(solutes[n])

    for i in range(len(nsolute)):
        for sol in tqdm.tqdm(range(nsolute[i])):
            solvent.place_solute_random(solutes[i])

    from pathlib import Path

    for p in Path(".").glob("step*"):
        p.unlink()
