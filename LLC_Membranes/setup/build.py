#! /usr/bin/env python3

import numpy as np
import argparse
from LLC_Membranes.setup.hexagonal_build import BuildHexagonal
from LLC_Membranes.setup.bcc_build import BicontinuousCubicBuild


def initialize():

    # TODO: make a yaml

    parser = argparse.ArgumentParser(description='Build LLC unit cell')

    parser.add_argument('-phase', '--phase', default='gyroid', type=str, help='Liquid crystal phase to build (HII, QI '
                                                                              'etc.)')
    parser.add_argument('-b', '--build_monomer', type=str, nargs='+', help='Name of single monomer'
                        'structure file (.gro format) used to build full system')
    parser.add_argument('-o', '--out', default='initial.gro', help='Name of output .gro file for full system')

    # System Geometry (phase agnostic)
    parser.add_argument('-r', '--pore_radius', default=.6, type=float, help='Initial Pore Radius (nm)')
    parser.add_argument('-seed', '--random_seed', default=False, type=int, help='Random seed for column shift. Set this'
                                                                                'to reproduce results')
    parser.add_argument('-box', '--box_lengths', nargs='+', type=float, help='Length of box vectors [x y z]')
    parser.add_argument('-angles', '--angles', nargs='+', default=[90, 90, 60], help='Angles between'
                        'box vectors')  # change this to None probably then set default angles based on phase

    # HII phase paramters
    parser.add_argument('-nc', '--ncolumns', default=5, type=int, help='Number of columns used to build each pore')
    parser.add_argument('-p', '--p2p', default=4.5, type=float, help='Initial pore-to-pore distance (nm)')
    parser.add_argument('-n', '--nopores', default=4, type=int, help='Number of pores (only works with 4 currently)')
    parser.add_argument('-d', '--dbwl', default=.37, type=float, help='Distance between vertically stacked monomers'
                                                                      '(nm)')
    parser.add_argument('-m', '--monomers_per_column', default=20, type=int, help='Number of monomers to stack in each'
                                                                                  'column')
    parser.add_argument('-t', '--tilt', default=0, type=float, help='Tilt angle of monomer with respect to xy plane (degrees)')
    parser.add_argument('-L', '--correlation_length', type=float, help='Length over which distance correlation between'
                                                                       'stacked monomers persists (nm)')
    parser.add_argument('--no_column_shift', action="store_false", help="Do not randomly shift columns")

    parser.add_argument('-Lvar', default=0.1, type=float, help='Variance in z position of monomer heads (nm)')
    parser.add_argument('-pd', '--parallel_displaced', default=0, type=float, help='Angle of wedge formed between line'
                        'extending from pore center to monomer and line from pore center to vertically adjacent monomer'
                                                                                   'head group.')
    parser.add_argument('-mf', '--mol_frac', nargs='+', default=[1.], type=float, help='If using the -random flag, this gives'
                        'the relative amount of each monomer to add. List the fractions in the same order as'
                        '-build_monomer')

    # Bicontinuous cubic structure
    parser.add_argument('-g', '--grid', default=50, type=int, help='Number of sections to break grid into when '
                                                                   'approximating the chosen implicit function')
    parser.add_argument('-dens', '--density', default=1.1, type=float, help='Density of system (g/cm3)')
    parser.add_argument('-c', '--curvature', default=1, type=float,
                        help='> 0 : QI phase (positive mean curvature), < 0, QII phase (negative mean'
                             'curvature). Determines whether the phase is normal or inverted')
    parser.add_argument('-wt', '--weight_percent', default=77.1, type=float,
                        help='Weight %% of monomer in membrane')
    parser.add_argument('-sol', '--solvent', default='glycerol', type=str,
                        help='Name of solvent mixed with monomer')
    parser.add_argument('-shift', '--shift', default=0, type=float,
                        help='Shift position of head group shift units '
                             'in the direction opposite of the normal vector at that point')
    parser.add_argument('-plot', '--plot_grid', action="store_true", help='Plot the grid of points used to defined the '
                                                                          'BCC surface')

    return parser


class PhaseError(Exception):
    """ Raised if invalid phase specified """

    def __init__(self, message):

        super().__init__(message)


if __name__ == "__main__":

    args = initialize().parse_args()

    # Choose phase and space group
    acceptable_hII_names = ['hii', 'h2', 'hexagonal', 'colh']  # lowercase so I can make input case insensitive
    acceptable_qI_names = ['gyroid', 'ia3d', 'schwarzd', 'pn3m', 'sphere']
    vague_names = ['q1', 'qi']  # phases by these names do not give enough information to determine the structure

    phase = args.phase.lower()

    if phase in acceptable_hII_names:
        phase = 'h2'
    elif phase in acceptable_qI_names:
        phase = 'q1'
    elif phase in vague_names:
        raise PhaseError('The phase you specified is not specfic enough to build an exact structure. For example, '
                         'QI does not distinguish between the gyroid and SchwarzD space groups.')
    else:
        raise PhaseError("'%s' does not specify a valid structure to build" % phase)

    if phase == 'h2':  # manipulate BuildHexagonal class (of build_hexagonal.py) with user-defined parameters

        if not args.build_monomer:
            build_monomer = 'NAcarb11V'  # default monomer for HII phase
        else:
            if type(args.build_monomer) is list:
                build_monomer = [i.split('.')[0] for i in args.build_monomer]
            else:
                build_monomer = args.build_monomer.split('.')[0]  # remove file extension if there is one

        correlation = False
        if args.correlation_length is not None:
            correlation = True

        if args.random_seed:
            np.random.seed(args.random_seed)
            seeds = np.random.randint(0, 4294967295, size=int(
                args.nopores * args.ncolumns))  # upper bound limit for numpy randint: see https://stackoverflow.com/questions/30721703/generate-random-integer-without-an-upper-bound

        system = BuildHexagonal(build_monomer, args.nopores, args.p2p, args.angles[2], args.pore_radius)

        system.reorient_monomer()  # align monomer with xy plane and orient along x-axis

        # system.align_plane()  # align monomer head group with xy plane
        # system.translate_to_origin()  # move monomer to origin for rotation
        # system.align_with_x()  # align vector from benzene ring to carboxylate with x axis

        wedge_theta = 360 / args.ncolumns  # rotation between laterally adjacent monomers (angle defining slice)
        for i in range(args.nopores):
            start_theta = 0
            thetas = [start_theta + x*wedge_theta for x in range(args.ncolumns)]
            for j in range(args.ncolumns):
                if args.random_seed:
                    np.random.seed(seeds[i * args.ncolumns + j])
                z = np.linspace(0, args.dbwl*args.monomers_per_column - args.dbwl, args.monomers_per_column)
                system.build_column(i, z, thetas[j], correlation=correlation, var=args.Lvar,
                                    correlation_length=args.correlation_length, pd=args.parallel_displaced,
                                    random_shift=args.no_column_shift, mole_fraction=args.mol_frac)

        system.reorder()
        print(system.LC[0].residues)
        exit()
        
        #system.exchange_ion()

        if args.box_lengths:
            a, b, c = args.box_lengths
        else:
            a, b, c = [2*args.p2p, 2*args.p2p, args.dbwl*args.monomers_per_column]  # logical choices for symmetry

        alpha, beta, gamma = [angle * (np.pi / 180) for angle in args.angles]

        V = a * b * c * np.sqrt(
            1 - np.cos(alpha) ** 2 - np.cos(beta) ** 2 - np.cos(gamma) ** 2 + 2 * np.cos(alpha) * np.cos(beta) * np.sin(
                gamma))  # volume of unitcell

        A = np.array([a, 0, 0])  # vector in x direction
        B = np.array([b*np.cos(gamma), b*np.sin(gamma), 0])  # vector in y direction
        C = np.array([c*np.cos(beta), c*((np.cos(alpha) - np.cos(gamma)*np.cos(beta))/np.sin(gamma)), V / (a*b*np.sin(gamma))])  # vector in z direction

        box = np.vstack((A, B, C))

        system.write_gro(args.out, box)

    elif phase == 'q1':

        if not args.build_monomer:
            build_monomer = 'Dibrpyr14'
        else:
            build_monomer = args.build_monomer.split('.')[0]  # get rid of file extension if there is one

        space_group = args.phase

        if not args.box_lengths:
            box_dimensions = 10.  # a default box size
        else:
            box_dimensions = args.box_lengths

        system = BicontinuousCubicBuild(build_monomer, space_group, box_dimensions, args.weight_percent, args.density)

        system.gen_grid(args.grid, args.curvature, plot=args.plot_grid)

        system.determine_monomer_placement(r=0.4)

        system.place_monomers(shift=args.shift)
        system.reorder()

        system.write_final_configuration()
