#!/usr/bin/env python

import argparse
import yaml
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from LLC_Membranes.analysis import hbonds, coordination_number


def initialize():

    parser = argparse.ArgumentParser(description='Model a continuous time random walk')

    # Pass yaml (preferred)
    parser.add_argument('-y', '--yaml', default=None, help='Name of configuration file. This is the preferred way to'
                                                           'pass parameters since it is easily reproducible. There are'
                                                           'also a lot of parameters')

    # args for everything in the yaml. These are ignored if yaml is passed
    # MD trajectory control
    parser.add_argument('-t', '--trajectory', default='PR_nojump.xtc', help='Path to input file.')
    parser.add_argument('-g', '--gro', default='em.gro', help='Name of .gro coordinate file.')
    parser.add_argument('-r', '--res', default='MET', help='Name of residue')

    # association parameters
    parser.add_argument('-rc', '--coordinated_residue', default=None, help='Name of residue coordinated to residue_')
    parser.add_argument('-ac', '--coordinated_atoms', default=None, nargs='+', help='Name of residue coordinate to '
                        'residue')
    parser.add_argument('-ca', '--catoms', default=None, nargs='+', help='Name of atoms to calculate coordination '
                        'number with respect to. The center of mass will be used')
    parser.add_argument('-ta', '--atype', default=None, help='Element name of atoms of which you want coordination '
                                                             'number')
    parser.add_argument('-tc', '--coordinated_type', default=None, help='Element name of coordinated atoms')
    parser.add_argument('-cut', default=0.25, type=float, help='Maximum distance between pairs where they are considered'
                                                              'coordinated (nm)')

    # hbond parameters

    parser.add_argument('-hr', '--hresidues', nargs='+', default=['HII'], help='Residues to include in h-bond search')
    parser.add_argument('-ha', '--hatoms', action='append', nargs='+', help='Atoms to include for each '
                        'residue. Each list of atoms must be passed with a separate -a flag for each residue in'
                        'args.residues')
    parser.add_argument('-hd', '--hdistance_cut', default=.35, type=float, help='Maximum distance between acceptor and'
                                                                                'donor atoms')
    parser.add_argument('-hangle', '--hangle_cut', default=30, type=float, help='Maximum DHA angle to be considered an '
                                                                                'H-bond')
    parser.add_argument('-hacc', '--acceptors', default=False, nargs='+', help='If you only want hbonds with acceptor'
                        'atoms of a residue, use this flag. 1 for this condition, otherwise 0. Input as a list'
                        ' of 0 and 1 in the same order as args.residues')
    parser.add_argument('-hdonors', '--donors', default=False, nargs='+', help='If you only want hbonds with donor'
                        'atoms of a residue, use this flag. 1 for this condition, otherwise 0. Input as a list'
                        ' of 0 and 1 in the same order as args.residues')

    return parser


class States:

    def __init__(self, traj, gro, res, **kwargs):
        """ Identify discrete states in a trajectory

        :param traj: name of GROMACS trajectory file (.xtc or .trr)
        :param gro: name of GROMACS coordinate file (.gro)
        :param res: name of residue whose states will be tracked
        :param kwargs: dictionaries of tuning paramters for states

        :type traj: str
        :type gro: str
        :type res: str
        :type kwargs: dict
        """

        # Read in all the parameters
        if not kwargs['association_params']:
            # default paramters for coordination_number.System. Can be modified with kwargs
            self.association_params = {'coordinated_residue': None, 'atoms': None, 'coordinated_atoms': None,
                                       'type': None, 'coordinated_type': None, 'begin': 0, 'end': -1, 'skip': 1,
                                       'com': True}
        else:
            self.association_params = kwargs['association_params']

        if not kwargs['hbond_params']:
            self.hbond_params = {'acceptors': False, 'residues': ['HII', 'MET'], 'donors': False, 'atoms': False,
                                 'angle_cut': 30, 'distance_cut': 0.35}
        else:
            self.hbond_params = kwargs['hbond_params']

        print("Loading trajectory...", end='', flush=True)
        self.t = md.load(traj, top=gro)
        print("Done!")

        print("Identifying hydrogen bonds...", end='', flush=True)
        self.hbonds = hbonds.System(traj, gro, t=self.t)
        self.identify_hydrogen_bonds()
        print("Done!")

        self.association = coordination_number.System(traj, gro, residue=res,
                                                      coordinated_residue=self.association_params['coordinated_residue'],
                                                      atoms=self.association_params['atoms'],
                                                      coordinated_atoms=self.association_params['coordinated_atoms'],
                                                      type=self.association_params['type'],
                                                      ctype=self.association_params['coordinated_type'],
                                                      com=self.association_params['com'], t=self.t)

    def identify_hydrogen_bonds(self):
        """ Use hbonds.System to identify hydrogen bonds between the residue and the specified set of atoms
        """

        residues = self.hbond_params['residues']
        if not self.hbond_params['acceptors']:
            acceptors = [False for _ in residues]
        else:
            acceptors = [bool(int(i)) for i in self.hbond_params['acceptors']]

        if not self.hbond_params['donors']:
            donors = [False for _ in residues]
        else:
            donors = [bool(int(i)) for i in self.hbond_params['donors']]

        if not self.hbond_params['atoms']:
            atoms = [['all'] for _ in residues]  # a default value
        else:
            atoms = self.hbond_params['atoms']

        while len(atoms) != len(residues):
            atoms.append(['all'])

        for i, r in enumerate(residues):
            self.hbonds.set_eligible(r, atoms[i], acceptor_only=acceptors[i], donors_only=donors[i])

        self.hbonds.identify_hbonds(self.hbond_params['distance_cut'], self.hbond_params['angle_cut'])


if __name__ == "__main__":

    args = initialize().parse_args()

    if not args.yaml:

        association_params = {'coordinated_residue': args.coordinated_residue, 'atoms': args.catoms,
                              'coordinated_atoms': args.coordinated_atoms, 'type': args.type,
                              'coordinated_type': args.coordinated_type,
                              'begin': 0, 'end': -1, 'skip': 1, 'com': True}

        hbond_params = {'acceptors': args.acceptors, 'residues': args.hresidues, 'donors': args.donors,
                        'atoms': args.hatoms, 'angle_cut': args.hangle_cut, 'distance_cut': args.hdistance_cut}

        traj, gro, res = args.trajectory, args.gro, args.res

    else:  # preferred

        with open(args.yaml, 'r') as yml:
            cfg = yaml.load(yml)

        trajectory = cfg['trajectory']
        traj, gro, res = trajectory['trajectory'], trajectory['gro'], trajectory['res']

        association_params = cfg['association_params']
        hbond_params = cfg['hbond_params']

    states = States(traj, gro, res, association_params=association_params, hbond_params=hbond_params)
