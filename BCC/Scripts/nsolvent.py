#!/usr/bin/env python

import argparse
import bcc_class

def initialize():

    parser = argparse.ArgumentParser(description='Create bicontinuous cubic membrane structure')  # allow input from user

    parser.add_argument('-b', '--build_mon', default='Dibrpyr14.gro', type=str, help='Name of monomer to build with')
    parser.add_argument('-d', '--dim', default=10, type=float, help='Unit cell dimension (length of x, y and z vector)')
    parser.add_argument('-dens', '--density', default=1.1, type=float, help='Density of system (g/cm3)')
    parser.add_argument('-wt', '--weight_percent', default=77.1, type=float, help='Weight %% of monomer in membrane')
    parser.add_argument('-sol', '--solvent', default='glycerol', type=str, help='Name of solvent mixed with monomer')

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = initialize()

    mw = {'glycerol': 92.09382}

    LC = bcc_class.LC(args.build_mon)  # get all properties of the liquid crystal

    # figure out how many solvent molecules are needed to achieve a certain weight percent
    NA = 6.022 * 10 ** 23  # avogadros number
    mass = LC.MW / NA  # mass of a single monomer (g)
    V = args.dim ** 3 * 10 ** -21  # volume of unit cell (cm ^ 3)

    mass = mw[args.solvent] / NA
    solv_mass = args.density * V * (100 - args.weight_percent) / 100
    nsol = int(solv_mass / mass)

    print(nsol)
