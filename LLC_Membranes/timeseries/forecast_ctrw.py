#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw
from LLC_Membranes.timeseries.ctrwsim import CTRW
from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
import sys
import yaml


def initialize():

    parser = argparse.ArgumentParser(description='Model a continuous time random walk')

    parser.add_argument('-y', '--yaml', help='Name of yaml configuration file with all of the following arguments'
                                             'already in it.')

    # MD trajectory control
    parser.add_argument('-t', '--trajectory', default='PR_nojump.xtc', help='Path to input file.')
    parser.add_argument('-g', '--gro', default='em.gro', help='Name of .gro coordinate file.')
    parser.add_argument('-b', '--begin', default=0, help='First trajectory frame used for analysis.')
    parser.add_argument('-e', '--end', default=-1, help='Last trajectory frame used for analysis.')
    parser.add_argument('-step', '--step', default=1, help='Include every "step" frames')
    parser.add_argument('-ma', '--moving_average', default=False, type=int, help='Calculate a moving average of the '
                        'center of mass coordinate')

    # loading saved objects (not part of yaml)
    parser.add_argument('-load', '--load', default=False, help='Specify name of pickled object to load')

    # restrict to pores parameters
    parser.add_argument('-restrict', '--restrict_to_pores', action='store_true', help='Only look at residues which'
                        'stay in the pore (based on last simulation frame)')
    parser.add_argument('--pore_defining_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Atoms'
                        'used to define pore centers')
    parser.add_argument('--pore_defining_residue', default='HII', type=str, help='residue to which pore_defining atoms '
                                                                                 'belong')
    parser.add_argument('-radius', '--pore_radius', default=1.48, type=float, help='Defined pore radius (nm)')

    parser.add_argument('-r', '--residue', default='ETH', type=str, help='Name of residue whose diffusivity we want')
    parser.add_argument('-atoms', nargs='+', help='Name of atoms whose collective diffusivity is desired')
    parser.add_argument('-nbins', default=25, type=int, help='Number of bins to bin hop and dwell distributions into')
    parser.add_argument('-bp', '--breakpoint_penalty', default=0.25, type=float, help='Cost function penalty when '
                                                                                      'determing break points')

    # ctrw simulation
    parser.add_argument('-ntsim', '--ntrajsim', default=1000, type=int, help='Number of trajectories to simulate')
    parser.add_argument('--update', action="store_true", help="update database with all parameters calculated")
    parser.add_argument('--ensemble', action="store_true", help="Calculate ensemble average MSD. If false, will do"
                                                                "time-averaged MSD")

    # bootstrapping
    parser.add_argument('-nboot', '--nboot', default=200, type=int, help='Number of bootstrap trials to be run')
    # parser.add_argument('-f', '--frontfrac', default=0.05, type=float, help='Where to start fitting line on msd curve')
    # parser.add_argument('-F', '--fracshow', default=.2, type=float, help='Percent of graph to show, also where to stop '
    #                     'fitting line during diffusivity calculation')
    # parser.add_argument('-a', '--axis', default='xyz', type=str, help='Which axis to compute msd along')
    # parser.add_argument('-nboot', default=200, type=int, help='Number of bootstrap trials for error estimation')
    # # | not implemented |
    # # v                 v
    # parser.add_argument('--restrict_to_pores', action="store_true", help='Only look at residue within pores of HII'
    #                                                                      'membrane')
    # parser.add_argument('-radius', default=1, type=float, help='Radius of pores. Anything greater than this distance'
    #                     'from the pore center will not be included in calculation')

    return parser


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.yaml:
        with open(args.yaml, 'r') as yml:
            cfg = yaml.load(yml)
    else:
        sys.exit('Argparse arguments (except --load and --yaml) are no longer supported. Please make a .yaml file.')

    if args.load:

        sys = file_rw.load_object(args.load)

    else:

        sys = SFBMParameters(cfg['trajectory'], cfg['gro'], cfg['residue'], start=cfg['begin'], end=cfg['end'],
                             step=cfg['step'], ma=cfg['ma'], nmodes=cfg['nmodes'])

        sys.calculate_solute_partition(spline=cfg['spline'], membrane_residue=cfg['membrane_residue'],
                                       r=cfg['pore_cut'])

        sys.hops_and_dwells(penalty=cfg['breakpoint_penalty'])

        sys.fit_distributions(nbins=cfg['bins'], nboot=cfg['nboot'], plot=True, show=False, save=True)

        sys.estimate_hurst()

        if args.update:
            sys.update_database()

        if sys.nmodes > 1:
            sys.determine_transition_matrix()

        sys.t = None  # save memory in pickled file this info is no longer needed. Can write function to reload
        sys.com = None
        file_rw.save_object(sys, 'forecast_%s_%dstate.pl' % (cfg['residue'], sys.nmodes))

    # simulate ntrajsim trajectories for same length as MD
    random_walks = CTRW(10000, args.ntrajsim, nmodes=sys.nmodes, dt=sys.dt, hop_dist='fbm', dwell_dist='power',
                        transition_matrix=sys.transition_matrix if sys.nmodes > 1 else None)

    random_walks.generate_trajectories(fixed_time=True, distributions=(sys.alpha_distribution,
                                       sys.hop_sigma_distribution, sys.hurst_distribution), discrete=True, ll=1)
    random_walks.calculate_msd(ensemble=args.ensemble)

    if args.ensemble:  # Ensemble-averaged MSD

        random_walks.bootstrap_msd(fit_power_law=True)
        random_walks.plot_msd(plot_power_law=True, show=False)
        # sys.update_database(type='msd', data=[1000 * (x / sys.time[-1]) for x in random_walks.final_msd])
        if args.update:
            sys.update_database(type='msd_ensemble', data=random_walks.final_msd)

    else:  # Time-averaged MSD

        random_walks.bootstrap_msd(fit_linear=False)
        random_walks.plot_msd(show=False, end_frame=8000)  # get data up to 400 ns
        if args.update:
            sys.update_database(type='msd_time_average', data=random_walks.final_msd)
