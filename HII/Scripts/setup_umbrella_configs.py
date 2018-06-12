#!/usr/bin/env python

import numpy as np
import place_solutes
import argparse
import genmdp
import os


def initialize():

    parser = argparse.ArgumentParser(description='Add specified amount of solvent to box')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Coordinate file to add solutes to')
    parser.add_argument('-n', '--n_configs', default=12, type=int, help='Number of configurations to generate')
    parser.add_argument('-d', '--spacing', type=float, default=0.2, help='Spacing between restraints')
    parser.add_argument('-s', '--solute', default='ETH', help='name of solute residue to add')
    parser.add_argument('-l', '--layers', default=20, type=int, help='number of layers in initial configuration')
    parser.add_argument('-p', '--pores', default=4, type=int, help='number of pores to add solutes to')
    parser.add_argument('-frac', default=0.5, type=float, help='Fraction into membrane to start placing solutes')
    parser.add_argument('-T', '--temp', default=300, type=float, help='Temperature to run simulations at')
    parser.add_argument('-o', '--output', default='umbrella')
    parser.add_argument('-r', '--ref', default=18, help='Layer number for reference com in each pore. The default value'
                        'of 18 will calculate the center of mass in each pore defined by the monomers that make up the'
                        '18th layer.')
    parser.add_argument('-L', '--sim_length', default=10000, type=int, help='Simulation length for each umbrella (ps)')
    parser.add_argument('-f', '--frames', default=100, type=int, help='Number of simulation frames to record. Will '
                        'determine nstxout, nstvout, and nstfout. Which means they will have the same value')
    parser.add_argument('-dt', default=0.002, type=float, help='Simulation time step (ps)')
    parser.add_argument('--mdp', action="store_true", help='Only write .mdp file. Do not create any configurations')
    parser.add_argument('-k', '--force_constant', type=float, help='Override force constant calculation and manually set force '
                                                       'constant (kJ / mol / nm^2)')

    args = parser.parse_args()

    return args


def create_pull_groups(system, solute, layer, nlayers, pores, ref_resname='HII', out='index.ndx', write=True):
    """
    :param system: object created by place_solutes.Solvent(). Will be used to make reference com groups
    :param solute: object created by place_solutes.Solute(). Will become a pull group
    :param layer: layer number which will be used as reference for com, starting at 1 (int)
    :param nlayers: number of layers per pore (int)
    :param pores: number of pores in system (int)
    :param resname: name of residue that will be restrained
    :param ref_resname: name of reference residue, some of whose atoms will define reference center of mass
    :param out: name of output index file
    :return: names of pull groups

    NOTE: currently this assumes one residue being pulled per pore. Each residue is a pull group and has its own
    center of mass pull group

    """

    ref_groups = []
    res_groups = []
    for i in range(pores):
        ref_groups.append('ref_%d' % (i + 1))
        res_groups.append('%s%d' % (solute.resname, i + 1))

    ref_res = [a.index for a in system.t.topology.atoms if a.residue.name == ref_resname]

    atoms_per_pore = len(ref_res) // pores
    atoms_per_layer = atoms_per_pore // nlayers

    if write:
        with open('%s' % out, 'a') as f:
            for i in range(pores):
                if i == 0:
                    f.write('[ %s ]\n' % ref_groups[i])
                else:
                    f.write('\n[ %s ]\n' % ref_groups[i])
                n = i * atoms_per_pore + (layer - 1)*atoms_per_layer
                for j in range(n, n + atoms_per_layer):
                    if (j - n) % 10 == 0 and (j - n) != 0:
                        f.write('%d\n' % (ref_res[j] + 1))  # GROMACS indexing starts at 1
                    else:
                        f.write('%d ' % (ref_res[j] + 1))
                f.write('\n')  # space between sections

            # written as a separate loop purely for readability in the index file
            for i in range(pores):
                f.write('\n[ %s ]\n' % res_groups[i])
                n = i * solute.natoms
                for j in range(n, n + solute.natoms):
                    if (j - n) % 10 == 0 and (j - n) != 0:
                        f.write('%d\n' % (system.positions.shape[0] + j + 1))  # GROMACS indexing starts at 1
                    else:
                        f.write('%d ' % (system.positions.shape[0] + j + 1))
                f.write('\n')

    return ref_groups, res_groups


if __name__ == "__main__":

    args = initialize()

    os.environ["GMX_MAXBACKUP"] = "-1"

    if args.force_constant:
        force_constant = args.force_constant
    else:
        kb = 1.381e-23 * 6.022e23 / 1000  # boltzmann constant kJ/(mol*K)
        beta = 1 / (kb * args.temp)  # mol / kJ
        force_constant = 1 / (beta * args.spacing ** 2)  # kJ / (mol*nm^2) -- force constant given a desired RMSD

    system = place_solutes.Solvent(args.gro)  # system to which to add solute
    solute = place_solutes.Solute(args.solute)  # object with all relevant solute properties

    # create index file with pull groups
    write = True
    if args.mdp:
        write = False
    ref_groups, res_groups = create_pull_groups(system, solute, args.ref, args.layers, args.pores, write=write)

    # create .mdp file
    nst = int((args.sim_length / args.dt) / args.frames)
    mdp = genmdp.SimulationMdp(args.gro, T=args.temp, length=args.sim_length, tau_p=1, nstxout=nst, nstvout=nst,
                               nstfout=nst)
    mdp.write_npt_mdp(out='pull')  # will write pull.mdp for an npt simulation (no pull parameters added yet)
    mdp.add_pull_groups(ref_groups, res_groups, force_constant, 0, 'pull.mdp')  # add pull parameters to pull.mdp

    if args.mdp:
        exit()

    # place solutes
    zbox = system.box_vectors[2, 2]
    z = np.linspace(zbox * args.frac, zbox*args.frac + args.n_configs*args.spacing - args.spacing, args.n_configs)

    for i, d in enumerate(z):
        if i != 0:  # save almost negligible amount of time..but why waste any?
            system = place_solutes.Solvent(args.gro)  # system to add solute to
        system.place_solute_pores(solute, d, layers=args.layers, pores=args.pores)
        # system.freeze_ndx(res=solute.resname)
        # system.energy_minimize(1000)  # energy minimize final configuration a bit longer
        system.write_config('%s_%s.gro' % (args.output, i))




