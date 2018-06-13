#! /usr/bin/env python

import mdtraj as md
import argparse
import Structure_char
import numpy as np
from llclib import physical
import matplotlib.pyplot as plt
import matplotlib as mpl
import os.path as path


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the number density of selected groups along axis')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='traj_whole.xtc', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop doing calculations')
    parser.add_argument('-bins', default=100, type=int, help='Number of bins to use')
    parser.add_argument('-m', '--multi', nargs='+', help='Overlay the density of each region with the results from '
                                                         'other trajectories')
    parser.add_argument('-s', '--solvate', action="store_true",
                        help='If the system is solvated, plot the number density of water as well')
    parser.add_argument('-l', '--load', help='Name of compressed .npz to load')

    args = parser.parse_args()

    return args


def duplicate(pos, box):
    """
    Duplicate a set of positions periodically once in the +/- xy directions
    :param pos: xyz positions of a set of coordinates to be duplicated periodically
    :param box: box vectors in mdtraj format (t.unitcell_vectors : [nframes, 3, 3]) for every frame
    :return: Periodically duplicated system
    """

    n = pos.shape[1]  # number atoms in original unit cell

    p = np.zeros([pos.shape[0], n*9, 3])  # will hold periodically duplicated system
    p[:, :n, :] = pos

    # x-direction
    for t in range(pos.shape[0]):
        p[t, n:2*n, :] = pos[t, :, :] + box[t, 0, :]
        p[t, 2*n:3*n, :] = pos[t, :, :] - box[t, 0, :]

    # y-direction
    n *= 3
    for t in range(pos.shape[0]):
        p[t, n:2*n, :] = p[t, :n, :] + box[t, 1, :]
        p[t, 2*n:3*n, :] = p[t, :n, :] - box[t, 1, :]

    return p


if __name__ == "__main__":
    
    args = initialize()  # parse the args

    # run : regional_density.py - m offset_disordered_regional_density.npz offset_regional_density.npz
    # layered_disordered_regional_density.npz layered_regional_density.npz solvated_regional_density.npz

    # mpl.style.use('seaborn')

    # regions = ['Tails', 'Head Groups', 'Sodium']
    regions = ['Tails', 'Head Groups', 'Sodium']

    if args.multi:

        # colors = ['blue', 'red', 'green', 'xkcd:orange']
        # colors = ['xkcd:red', 'xkcd:green', 'blue', 'xkcd:yellow']
        colors = ['xkcd:blue', 'xkcd:olive', 'xkcd:orangered', 'xkcd:magenta', 'xkcd:gold']
        names = ['Ordered Parallel Displaced', 'Ordered Sandwiched', 'Disordered Sandwiched',
                 'Disordered Parallel Displaced', 'Solvated Parallel Displaced']
        # names = ['Dry', 'Solvated']
        # colors = ['xkcd:orange', 'xkcd:blue', 'xkcd:orange']

        n = len(args.multi)
        system = np.load(args.multi[0])

        # It is assumed that all of the data uses the same number of bins with the same bin width
        r = system['r']
        bin_width = system['bw']

        results = np.zeros([n, len(regions), len(r)])
        # results = np.zeros([n, 4, len(r)])

        # results[0, :3, :] = system['results']

        for i in range(n - 1):
            system = np.load(args.multi[i])
            results[i, :, :] = system['results']

        system = np.load(args.multi[-1])
        results[-1, :, :] = system['results'][:3, :]

        outline = np.zeros([4, n, r.shape[0]*2 + 2, 2])
        for i in range(len(regions)):
            plt.figure(i)
            for j in range(n):
                #plt.bar(r, results[j, i, :], bin_width, color=colors[j], alpha=1, label=names[j])

                half_width = bin_width / 2
                # outline = np.zeros([r.shape[0]*2 + 2, 2])
                outline[i, j, 0, 0] = r[0] - half_width
                outline[i, j, -1, 0] = r[-1] + half_width
                outline[i, j, 1:-1:2, 0] = r - half_width
                outline[i, j, 2:-1:2, 0] = r + half_width
                outline[i, j, 1:-1:2, 1] = results[j, i, :]
                outline[i, j, 2:-1:2, 1] = results[j, i, :]
                if i == 0:
                    plt.plot(outline[i, j, :-1, 0], outline[i, j, :-1, 1], color=colors[j], linewidth=2, label=names[j])
                else:
                    plt.plot(outline[i, j, 1:, 0], outline[i, j, 1:, 1], color=colors[j], linewidth=2,
                             label=names[j])
            # plt.title(regions[i], fontsize=14)
            plt.legend(fontsize=11)
            plt.ylabel('Component Number Density (number/nm$^2$)', fontsize=14)
            plt.xlabel('Distance from pore center, r (nm)', fontsize=14)
            plt.axes().tick_params(labelsize=14)
            # plt.ylim([0, 0.6])
            plt.tight_layout()
            plt.savefig("%s_density.png" % regions[i])

        plt.show()

        exit()

    if args.solvate:
        colors = ['xkcd:red', 'xkcd:green', 'blue', 'xkcd:yellow']
    else:
        colors = ['red', 'green', 'blue']

    if not args.load:

        print('Loading trajectory...', end="")
        t = md.load('%s' % args.traj, top='%s' % args.gro)[args.begin:args.end]
        print('done')

        box = t.unitcell_vectors
        nT = t.n_frames
        npores = 4
        r_max = 0

        if args.solvate:
            results = np.zeros([len(regions) + 1, args.bins])
        else:
            results = np.zeros([len(regions), args.bins])

        #keep = [a.index for a in t.topology.atoms if a.residue.name != 'HOH']  # everything kept if system not solvated
        components = ['C', 'C1', 'C2', 'C3', 'C4', 'C5']
        comp = [a.index for a in t.topology.atoms if a.name in components]

        p_centers = physical.avg_pore_loc(npores, t.xyz[:, comp, :])

        for i, reg in enumerate(regions):

            print('Calculating number density of %s region' % reg)

            pos = Structure_char.restrict_atoms(t, reg)  # restrict trajectory to region

            p = duplicate(pos, t.unitcell_vectors)  # duplicate things periodically

            equil = 0
            density, r, bin_width = Structure_char.compdensity(pos, p_centers, equil, box, pores=npores, buffer=0, nbins=args.bins)

            results[i, :] = density

            plt.bar(r, density, bin_width, color=colors[i], alpha=0.6, label=reg)

        np.savez_compressed("regional_density", results=results, r=r, bw=bin_width, box=t.unitcell_vectors)
        print('Arrays saved as density.npz')

        if args.solvate:

            print('Calculating number density of solvent')
            keep = [a.index for a in t.topology.atoms if a.residue.name == 'HOH' and a.name == 'O']
            pos = t.xyz[:, keep, :]
            p = duplicate(pos, t.unitcell_vectors)
            equil = 0
            density, r, bin_width = Structure_char.compdensity(pos, p_centers, equil, t.unitcell_vectors, pores=npores, buffer=0, nbins=args.bins)
            results[-1, :] = density
            plt.bar(r, density, bin_width, color='xkcd:gold', alpha=0.75, label='Water')

    else:

        d = np.load(args.load)

        bw = d['bw']
        r = d['r']
        results = d['results']
        box = d['box']

        if args.solvate:
            regions.append('Water')

        for i in range(results.shape[0]):
            plt.bar(r, results[i, :], bw, color=colors[i], alpha=0.75, label=regions[i])

    un_normalize = np.zeros_like(results)

    annulus_area = []
    dr = r[-1] - r[-2]  # width of annulus
    # r from above is calculated as the average of the edges of the bins. So I need to reverse that here
    for i in r:
        inner_r = i - (dr/2)
        outer_r = i + (dr/2)
        annulus_area.append(np.pi*(outer_r**2 - inner_r**2))

    for i in range(results.shape[0]):
        un_normalize[i, :] = results[i, :] * annulus_area

    r_cut = 0.6
    cut_index = 0
    while r[cut_index] < r_cut:
        cut_index += 1

    print('Percentage of sodium within %s nm of pore center: %2.2f %%' % (r_cut, 100*(sum(un_normalize[2, :cut_index])/sum(un_normalize[2, :]))))

    r_cut = 0.6
    cut_index = len(r) - 1
    while r[cut_index] > r_cut:
        cut_index -= 1

    print('Percentage of tails within %s nm of pore center: %2.2f %%' % (r_cut, 100*(sum(un_normalize[0, :cut_index])/sum(un_normalize[0, :]))))

    r_cut_pore = r_cut
    pore_cut_index = cut_index

    r_cut = 1
    cut_index = len(r) - 1
    while r[cut_index] > r_cut:
        cut_index -= 1

    print('Percentage of tails between %s and %s nm of pore center: %2.2f %%' % (r_cut_pore, r_cut,
        100*(sum(un_normalize[0, pore_cut_index:cut_index])/sum(un_normalize[0, :]))))

    print('Maximum of Head group region: r = %s' % r[np.argmax(results[1, :])])

    plt.legend(prop={'size':14}, loc=1)
    plt.ylabel('Component number density (count/nm$^3$)', fontsize=14)
    plt.xlabel('Distance from pore center (nm)', fontsize=14)
    plt.axes().tick_params(labelsize=14)
    plt.xlim(0, 1.5)
    plt.tight_layout()
    plt.savefig("regional_density.png")
    plt.show()
