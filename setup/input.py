#! /usr/bin/env python

"""
Write input files. Default settings are not included. Add an argument if you need to change a default that isn't shown
"""
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
import os
import argparse
from LLC_Membranes.setup import lc_class, bcc_class
from LLC_Membranes.setup.gentop import SystemTopology
from LLC_Membranes.setup.genmdp import SimulationMdp


def initialize():

    parser = argparse.ArgumentParser(description='Write .mdp files and topology files for various types of simulation')

    parser.add_argument('-T', '--title', default='Generic Molecular Dynamics Simulation', type=str, help='Simulation Title')
    parser.add_argument('-b', '--build_mon', default='NAcarb11V', type=str, help='Monomer structure used for build')
    parser.add_argument('-s', '--em_steps', default=5000, type=int, help='Steps to take during energy minimization')
    parser.add_argument('-e', '--ensemble', default='npt', type=str, help='Thermodynamic ensemble to put system in')
    parser.add_argument('-d', '--dt', default=0.002, type=float, help='time step (ps)')
    parser.add_argument('-l', '--length', default=1000, type=float, help='simulation length (ps)')
    parser.add_argument('-f', '--frames', default=500, type=int, help='number of frames')
    parser.add_argument('-p', '--pcoupltype', default='semiisotropic', type=str, help='Pressure Couple Type')
    parser.add_argument('--restraints', nargs='+', help='Specify residue for which restraints are specified')
    parser.add_argument('-x', '--xlink', help='Turn this to "on" if the the system is crosslinked', action="store_true")
    parser.add_argument('-c', '--coord', default='initial.gro', type=str, help='coordinate file of system to be simulated')
    parser.add_argument('-S', '--solvate', help='Specify this if the system has water so an extra line can be added to the '
                                                'topology', action="store_true")
    parser.add_argument('--temp', default=300, help='Specify temperature at which to run simulation')
    parser.add_argument('--mdp', action="store_true", help='Only the .mdp will be written if this option is specified')
    parser.add_argument('--barostat', default='berendsen', type=str, help='pressure coupling scheme to use')
    parser.add_argument('--genvel', default='yes', type=str, help='generate velocities according to a maxwell'
                                                                     'distribution')
    parser.add_argument('--bcc', action="store_true", help='Generate input files using bicontinuous cubic files')  # probably best to reorganize the repository
    parser.add_argument('--solvent', default='water', help='Name of solvent')
    parser.add_argument('--tau_t', default=0.1, type=float, help='Temperature coupling time constant')
    parser.add_argument('--tau_p', default=20, type=float, help='Pressure coupling time constant')
    parser.add_argument('-nx', '--nstxout', type=int, help='Frequency to output coordinates to trajectory file')
    parser.add_argument('-nv', '--nstvout', type=int, help='Frequency to output velocities to trajectory file')
    parser.add_argument('-nf', '--nstfout', type=int, help='Frequency to output forces to trajectory file')
    parser.add_argument('-ne', '--nstenergy', type=int, help='Frequency to output energy to energy file')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = initialize()

    # get output frequencies (important for controlling size of trajectory)
    if args.nstxout:
        nstxout = args.nstxout
    else:
        nstxout = int(args.length / (args.dt * args.frames))

    if args.nstvout:
        nstvout = args.nstvout
    else:
        nstvout = int(args.length / (args.dt * args.frames))

    if args.nstfout:
        nstfout = args.nstfout
    else:
        nstfout = int(args.length / (args.dt * args.frames))

    if args.nstenergy:
        nstenergy = args.nstenergy
    else:
        nstenergy = int(args.length / (args.dt * args.frames))  # output frequency

    mdp = SimulationMdp(args.coord, title=args.title, T=args.temp, em_steps=args.em_steps,
                        time_step=args.dt, length=args.length, p_coupling=args.pcoupltype,
                        barostat=args.barostat, genvel=args.genvel, restraints=args.restraints, xlink=args.xlink,
                        bcc=args.bcc, tau_p=args.tau_p, tau_t=args.tau_t, nstxout=nstxout, nstvout=nstvout,
                        nstfout=nstfout, nstenergy=nstenergy)

    mdp.write_em_mdp()  # write energy minimization .mdp without asking

    if args.ensemble == 'npt':
        mdp.write_npt_mdp()
    elif args.ensemble == 'nvt':
        mdp.write_nvt_mdp()

    top = SystemTopology(args.coord, restraints=args.restraints, xlink=args.xlink)
    top.write_top()
