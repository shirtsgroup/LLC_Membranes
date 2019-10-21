#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw
from LLC_Membranes.timeseries.ctrwsim import CTRW
from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
import sys
import yaml
import numpy as np


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

    parser.add_argument('-suffix', '--save_suffix', default=None, help='If not None, add this suffix to save name')

    return parser


def calculate_moving_average(series, n):
    """ Calculate moving average of a time series

    :param n: Number of previous points to average

    :type n: int
    """

    ret = np.cumsum(series, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]

    return ret[n - 1:] / n


if __name__ == "__main__":

    args = initialize().parse_args()

    if args.yaml:
        with open(args.yaml, 'r') as yml:
            cfg = yaml.load(yml)

        fit = cfg['fit']
        project = cfg['project']

    else:
        sys.exit('Argparse arguments (except --load and --yaml) are no longer supported. Please make a .yaml file.')

    if args.load:

        sys = file_rw.load_object(args.load)

    else:

        sys = SFBMParameters(fit['trajectory'], fit['gro'], fit['residue'], start=fit['begin'], end=fit['end'],
                             step=fit['step'], ma=fit['ma'], nmodes=fit['nmodes'])

        sys.calculate_solute_partition(spline=fit['spline'], membrane_residue=fit['membrane_residue'],
                                       r=fit['pore_cut'])

        # for x in [0.5, 0.75, 1.0, 1.25, 1.50]:
        #     sys.calculate_solute_partition(spline=fit['spline'], membrane_residue=fit['membrane_residue'],
        #                                    r=x)
        #     plt.plot(sys.time[:-24], calculate_moving_average(sys.partition.sum(axis=1), 25) / 24,
        #              label='r = %.2f nm' % x, lw=2)
        # plt.ylabel('Fraction of solutes in pore', fontsize=14)
        # plt.xlabel('Time (ns)', fontsize=14)
        # plt.legend()
        # plt.show()
        # exit()

        sys.hops_and_dwells(penalty=fit['breakpoint_penalty'])

        sys.fit_distributions(nbins=fit['bins'], nboot=fit['nboot'], plot=False, show=False, save=True,
                              dwell_distribution=fit['dwell_distribution'], hop_distribution=fit['hop_distribution'])

        sys.estimate_hurst(modes=fit['hurst_modes'], show=False)

        if args.update:
            sys.update_database()

        # if sys.nmodes > 1:
        #     sys.determine_transition_matrix()

        sys.t = None  # save memory in pickled file this info is no longer needed. Can write function to reload
        sys.com = None

        savename = 'forecast_%s_%dstate' % (fit['residue'], sys.nmodes)
        if args.save_suffix is not None:
            savename += args.save_suffix
        savename += '.pl'

        file_rw.save_object(sys, savename)

    # sys.fit_distributions(nbins=fit['bins'], nboot=10, plot=True, show=True, save=True,
    #                       dwell_distribution=fit['dwell_distribution'], hop_distribution=fit['hop_distribution'])

    # sys.estimate_hurst(modes=fit['hurst_modes'])
    # exit()
    exit()
    print('Generating SFBM realizations...')
    # simulate ntrajsim trajectories for same length as MD

    nsteps = int(project['length'] / sys.dt)
    random_walks = CTRW(nsteps, project['ntrajsim'], nmodes=sys.nmodes, dt=sys.dt, hop_dist=project['hop_dist'],
                        dwell_dist=project['dwell_dist'],
                        transition_count_matrix=sys.count_matrix if sys.nmodes > 1 else None)

    random_walks.generate_trajectories(fixed_time=True, distributions=(sys.dwell_parameters,
                                       sys.hop_parameters, sys.hurst_distribution), discrete=True,
                                       ll=sys.dwell_lower_limit, max_hop=sys.max_hop)  # change lower limit

    if project['ensemble']:
        print('Calculating ensemble averaged mean squared displacement...')
    else:
        print('Calculating time averaged mean squared displacement...')

    random_walks.calculate_msd(ensemble=project['ensemble'])

    if project['ensemble']:  # Ensemble-averaged MSD

        random_walks.bootstrap_msd(fit_power_law=True)
        random_walks.plot_msd(plot_power_law=True, show=False)
        # sys.update_database(type='msd', data=[1000 * (x / sys.time[-1]) for x in random_walks.final_msd])
        if args.update:
            sys.update_database(type='msd_ensemble', data=random_walks.final_msd)

    else:  # Time-averaged MSD

        random_walks.bootstrap_msd(fit_linear=False)
        random_walks.plot_msd(show=True, end_frame=int(project['frac_MSD_show'] * project['padding'] * nsteps))
        if args.update:
            sys.update_database(type='msd_time_average', data=random_walks.final_msd)
