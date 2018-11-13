#!/usr/bin/env python

import numpy as np
from LLC_Membranes.setup import place_solutes, genmdp
import argparse
import os
import subprocess


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
    parser.add_argument('-rl', '--ref_layer', default=3, help='Layer number for reference com in each pore. The default value'
                        'of 18 will calculate the center of mass in each pore defined by the monomers that make up the'
                        '18th layer.')
    parser.add_argument('-L', '--sim_length', default=10000, type=int, help='Simulation length for each umbrella (ps)')
    parser.add_argument('-f', '--frames', default=100, type=int, help='Number of simulation frames to record. Will '
                        'determine nstxout, nstvout, and nstfout. Which means they will have the same value')
    parser.add_argument('-dt', default=0.002, type=float, help='Simulation time step (ps)')
    parser.add_argument('--mdp', action="store_true", help='Only write .mdp file. Do not create any configurations')
    parser.add_argument('-k', '--force_constant', type=float, help='Override force constant calculation and manually set force '
                                                       'constant (kJ / mol / nm^2)')
    parser.add_argument('-r', '--ref', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], nargs='+', help='Atoms to use for'
                        'center of mass reference groups. Also dictates placement of solutes')

    args = parser.parse_args()

    return args


def calculate_force_constant(rmsd):
    """
    Calculate force constant require to hold a molecule in place with a desired root mean square deviation from the
    reference position
    :param rmsd: deviation from reference position (float, nm)
    """

    kb = 1.381e-23 * 6.022e23 / 1000  # boltzmann constant kJ/(mol*K)
    beta = 1 / (kb * args.temp)  # mol / kJ

    return 1 / (beta * rmsd ** 2)  # kJ / (mol*nm^2) -- force constant given a desired RMSD


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


def create_individual_pull_groups(system, solute, n_configs, layer, nlayers, pores,
                                  ref=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], out='index.ndx', write=True):
    """
    Give each molecule that is being pulled its own pull group. This assumes that the first molecule is placed at the
    center of mass (COM) of the first layer (defined by 'layer' parameter). The second molecule is place vertically
    above at the COM of the first and second layers. The third is at the COM of the second layer and the alternation
    continues from there.

    Naming convention for pull groups:
        - reference groups: [ ref_config-pore ] (for example, ref_1-2 refers to the reference group for the solute in
        the first configuration, in pore 2)
        - pull groups: [ solutename_config-pore ] (for example, ETH_2-4 referes to an ethanol molecule in the second
        configuration, in pore 4)

    :param system: object created by place_solutes.Solvent(). Will be used to make reference com groups
    :param solute: object created by place_solutes.Solute(). Will become a pull group
    :param n_configs : number of configurations
    :param layer: layer number which will be used as reference for 1st com, starting at 1 (int). Assumes the next ref
    group goes directly above.
    :param nlayers: number of layers per pore (int)
    :param pores: number of pores in system (int)
    :param resname: name of residue that will be restrained
    :param ref_resname: name of reference residue, some of whose atoms will define reference center of mass
    :param out: name of output index file
    :param write: if True, index file will be written (bool)
    :return: names of pull groups

    """

    ref_groups = []
    res_groups = []
    for j in range(n_configs):
        ref_groups.append([])
        res_groups.append([])
        for i in range(pores):
            ref_groups[j].append('ref_%d-%d' % (j + 1, i + 1))
            res_groups[j].append('%s_%d-%d' % (solute.resname, j + 1, i + 1))

    ref_res = [a.index for a in system.t.topology.atoms if a.name in ref]

    atoms_per_pore = len(ref_res) // pores
    atoms_per_layer = atoms_per_pore // nlayers

    if write:
        with open('%s' % out, 'a') as f:
            for j in range(n_configs):
                layers = 1  # number of layers included in COM calculation
                if j % 2 == 1:  # handles case when solute is placed between layers (uses COM of both layers)
                    layers = 2
                for i in range(pores):
                    f.write('\n[ %s ]\n' % ref_groups[j][i])
                    n = i * atoms_per_pore + (layer + j//2 - 1)*atoms_per_layer
                    for k in range(n, n + layers*atoms_per_layer):
                        if (k - n) % 10 == 0 and (k - n) != 0:
                            f.write('%d\n' % (ref_res[k] + 1))  # GROMACS indexing starts at 1
                        else:
                            f.write('%d ' % (ref_res[k] + 1))
                    f.write('\n')  # space between sections

            for j in range(n_configs):
                for i in range(pores):
                    f.write('\n[ %s ]\n' % res_groups[j][i])
                    n = i * solute.natoms
                    for k in range(n, n + solute.natoms):
                        if (k - n) % 10 == 0 and (k - n) != 0:
                            f.write('%d\n' % (system.positions.shape[0] + k + 1))  # GROMACS indexing starts at 1
                        else:
                            f.write('%d ' % (system.positions.shape[0] + k + 1))
                    f.write('\n')

    return ref_groups, res_groups


if __name__ == "__main__":

    args = initialize()

    os.environ["GMX_MAXBACKUP"] = "-1"

    if args.force_constant:
        force_constant = args.force_constant
    else:
        force_constant = calculate_force_constant(args.spacing)

    print('Initializing system...', flush=True, end='')
    system = place_solutes.Solvent(args.gro)  # system to which to add solute
    solute = place_solutes.Solute(args.solute)  # object with all relevant solute properties
    print('Success!')

    # find center of mass. Assumes all atoms in ref are the same. trace_pores will need to be updated to properly
    # do center of masses of different atoms
    print('Determing placement...', flush=True, end='')
    ref = [a.index for a in system.t.topology.atoms if a.name in args.ref]
    pore_spline = place_solutes.trace_pores(system.positions[ref, :], system.box_vectors[:2, :2], args.layers)

    # only need z placement values. place_solutes_pores will take care of xy position.
    z = np.zeros([args.pores, args.n_configs])
    n_layers = args.n_configs // 2 + 1  # number of layers from pore spline needed
    for i in range(args.pores):
        ndx = i * args.layers + args.ref_layer - 1
        layer_locations = pore_spline[ndx:(ndx + n_layers), 2]
        between_layers = np.array([(layer_locations[i] + layer_locations[i - 1])/2 for i in range(1, len(layer_locations))])
        z[i, ::2] = layer_locations[:z[0, ::2].size]
        z[i, 1::2] = between_layers
    print('Success!')
    print(z)
    exit()
    with open('centers.txt', 'w') as f:
        f.write('  pore 1 |  pore 2 |  pore 3 |  pore 4\n')
        for k in range(args.n_configs):
            for i in range(args.pores):
                f.write('  %.3f   ' % z[i, k])
            f.write('\n')
    print('Absolute starting center of mass recorded to centers.txt')

    # place solutes at equally spaced locations
    # zbox = system.box_vectors[2, 2]
    # z = np.linspace(zbox * args.frac, zbox*args.frac + args.n_configs*args.spacing - args.spacing, args.n_configs)

    print('Placing solutes...')
    if not args.mdp:
        for i, d in enumerate(z.T):
            if i != 0:  # save almost negligible amount of time..but why waste any?
                system = place_solutes.Solvent(args.gro)  # system to add solute to
            system.place_solute_pores(solute, d, layers=args.layers, pores=args.pores)
            system.freeze_ndx(res=solute.resname)
            system.energy_minimize(1000, freeze=True)  # energy minimize final configuration a bit longer
            system.write_config('%s_%s.gro' % (args.output, i + 1))
    print('Success!')

    print('Creating index file...', flush=True, end='')
    # create index file with pull groups
    write = True
    # if args.mdp:
    #     write = False

    # write default index groups using GROMACS (creates index.ndx by default)
    ps = subprocess.Popen(['echo', 'q'], stdout=subprocess.PIPE)
    ps2 = subprocess.Popen(['gmx', 'make_ndx', '-f', '%s_1.gro' % args.output], stdin=ps.stdout,
                           stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
    ps2.wait()

    ref_groups, res_groups = create_individual_pull_groups(system, solute, args.n_configs, args.ref_layer, args.layers,
                                                           args.pores, ref=args.ref, write=write)
    print('Success!')

    print('Creating .mdp files and atomic level input (.tpr) files...', flush=True, end='')
    nst = int((args.sim_length / args.dt) / args.frames)
    mdp = genmdp.SimulationMdp(args.gro, T=args.temp, length=args.sim_length, tau_p=1, nstxout=nst, nstvout=nst,
                               nstfout=nst)

    # write separate .mdp files for each system since pull groups are different, then grompp them
    for i in range(args.n_configs):
        mdp.write_npt_mdp(out='pull_%d' % (i + 1))  # write pull.mdp for an NPT simulation without pull parameters
        mdp.add_pull_groups(ref_groups[i], res_groups[i], force_constant, 0, 'pull_%d.mdp' % (i + 1))  # add pull params
        p = subprocess.Popen(['gmx', 'grompp', '-f', 'pull_%d.mdp' % (i + 1), '-p', 'topol.top', '-n', 'index.ndx',
                              '-c', '%s_%d.gro' % (args.output, i + 1), '-o', '%s_%d' % (args.output, i + 1)],
                             stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)
        p.wait()
    print('Success!')
