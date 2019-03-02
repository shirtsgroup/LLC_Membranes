#! /usr/bin/env python

import mdtraj as md
import argparse
import numpy as np
from LLC_Membranes.llclib import physical, topology
from LLC_Membranes.analysis import p2p
import matplotlib.pyplot as plt
import tqdm
import matplotlib as mpl
import os.path as path


def initialize():

    parser = argparse.ArgumentParser(description='Calculate the number density of selected groups along axis')

    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='traj_whole.xtc', help='Trajectory file (.trr, .xtc should work)')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame to begin calculations')
    parser.add_argument('-c', '--components', nargs='+', help="Region(s) or atoms to include in calculation. Special"
                                                              "groups include 'head groups' and 'tails'")
    parser.add_argument('-m', '--monomer', default='NaGA3C11', help='Name of monomer (no file extension)')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Frame to stop doing calculations')
    parser.add_argument('-bins', default=100, type=int, help='Number of bins to use')
    parser.add_argument('-s', '--solvate', action="store_true",
                        help='If the system is solvated, plot the number density of water as well')
    parser.add_argument('-l', '--load', help='Name of compressed .npz to load')

    args = parser.parse_args()

    return args


def grps(name):
    """ Return names of atoms making up specialized groups. This is specifically for NaGA3C11.

    TODO: incorporate these groups as annotations (maybe)

    :param name: specialized group names. Valid options are 'head groups' and 'tails'

    :type name: str

    :return: atom names that constitute specialized groups
    """

    if name == 'head groups':

        return ['C', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O3', 'O4']

    elif name == 'tails':

        return ['C7', 'C8', 'C9', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21',
                'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35',
                'C36', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42', 'C43', 'C44', 'C45', 'C46', 'C47', 'C48']


class Region(object):

    def __init__(self, traj, gro, components, monomer, begin=0, end=-1):
        """ Determine the density of components. Similar to a radial density, but instead of center of mass, the radial
        density of each atom in the region is averaged together.

        :param traj: name of trajectory to analyze (.xtc or .trr)
        :param gro: name of GROMACS coordinate file with same topology as traj
        :param component: list of atoms or name of groups of atoms
        :param begin: first frame to analyze
        :param end: last grame to analyze

        :type traj: str
        :type gro: str
        :type component: list
        :type begin: int
        :type end: int
        """

        special_groups = ['head groups', 'tails']

        if components[0].lower() in special_groups:
            self.components = grps(components[0].lower())
        else:
            self.components = components

        print('Loading trajectory...', end="")
        self.t = md.load('%s' % traj, top='%s' % gro)[begin:end]
        print('done')

        self.monomer = topology.LC('%s.gro' % monomer)

        self.r = None
        self.density = None

    def component_density(self, bins=50, spline=True, progress=True, npts_spline=10):

        self.r = np.zeros([bins])
        self.density = np.zeros([self.t.n_frames, bins])

        pore_defining_atoms = [a.index for a in self.t.topology.atoms if a.name in self.monomer.pore_defining_atoms
                               and a.residue.name in self.monomer.residues]

        pore_centers = physical.avg_pore_loc(4, self.t.xyz[:, pore_defining_atoms, :], self.t.unitcell_vectors,
                                             spline=spline, progress=progress, npts=npts_spline)

        print('Calculating component density')
        self.r, self.density = physical.compdensity(self.com, pore_centers, self.t.unitcell_vectors,
                                                    nbins=bins, spline=spline, cut=cut)


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


def compdensity(component, pore_centers, start, box, cut=1.5, pores=4, nbins=50, rmax=3.5, buffer=0.0):
    """ Measure the density of a component as a function of the distance from the pore centers

    :param component: the coordinates of the component(s) which you want a radial distribution of at each frame
    :param pore_centers: a numpy array of the locations of each pore center at each trajectory frame
    :param start: the frame number at which to start calculations (should be after equilibration)
    :param cut: cutoff distance for distance calculations. Will not count anything further than cut from the pore center
    :param pores: number of pores (int) default=4
    :param rmax: maximum distance from pore center to calculate density for, default = 3.5 nm
    :param buffer: percentage used to define the location of z planes between which component density will be computed,
           float, default = 0 (i.e. no buffer). Should be between 0 and 1. e.g. for 1 percent, use 0.01 as the buffer

    :type component: numpy.ndarray
    :type pore_centers: numpy.ndarray
    :type start: int
    :type cut: float
    :type pores: int
    :type rmax: float
    :type buffer: float

    :return: the density of "component" as a function the distance from the pore center. Also
             returns the calculated bin width for plotting
    """

    # Extract basic system information. It's important to follow the format of the component array to get it right

    tot_atoms = component.shape[1]  # the total number of components in a single frame
    n_ppore = tot_atoms // pores  # the total number of components in each pore

    nT = component.shape[0]
    zbox = np.mean(box[:, 2, 2])
    # Find the approximate max and minimum z values of the components based on the last frame

    zmax = np.max(component[-1, :, 2])
    zmin = np.min(component[-1, :, 2])
    thickness = zmax - zmin  # approximate membrane thickness

    # now find the maximum and minimum permissible z dimensions based on the buffer
    zmax_buff = zmax - thickness * buffer  # Could use buffer/2 depending on how you interpret what buffer % means
    zmin_buff = zmin + thickness * buffer

    density = np.zeros([nbins])  # number / nm^3
    for t in tqdm.tqdm(range(start, nT)):
        for p in range(pores):
            # narrow down the positions to those that are within 'cut' of at least one pore
            distances = np.linalg.norm(component[t, :, :2] - pore_centers[t, p, :], axis=1)
            d_sorted = np.sort(distances)
            # find where the distances exceed the cutoff
            stop = 0
            while d_sorted[stop] < cut:
                stop += 1

            hist, bin_edges = np.histogram(d_sorted[:stop], bins=nbins, range=(0, cut))  # the range option is necessary
            #  to make sure we have equal sized bins on every iteration

            density += hist

    density /= zbox * (nT - start)  # take average
    bin_width = cut / nbins

    # normalize based on area of anulus where bin is located
    r = np.zeros([nbins])
    normalization = []
    for i in range(nbins):
        normalization.append(np.pi * (bin_edges[i + 1] ** 2 - bin_edges[i] ** 2))
        density[i] /= np.pi * (bin_edges[i + 1] ** 2 - bin_edges[i] ** 2)
        r[i] = (bin_edges[i + 1] + bin_edges[i]) / 2

    return density, r, bin_width


if __name__ == "__main__":
    
    args = initialize()  # parse the args

    # run : regional_density.py - m offset_disordered_regional_density.npz offset_regional_density.npz
    # layered_disordered_regional_density.npz layered_regional_density.npz solvated_regional_density.npz

    # mpl.style.use('seaborn')

    # regions = ['Tails', 'Head Groups', 'Sodium']
    regions = ['Tails', 'Head Groups', 'Sodium']

    # if args.multi:
    if False:

        # colors = ['blue', 'red', 'green', 'xkcd:orange']
        # colors = ['xkcd:red', 'xkcd:green', 'blue', 'xkcd:yellow']
        colors = ['xkcd:blue', 'xkcd:olive', 'xkcd:orangered', 'xkcd:magenta', 'xkcd:gold']
        names = ['Ordered Parallel Displaced', 'Ordered Sandwiched', 'Disordered Sandwiched',
                 'Disordered Parallel Displaced', 'Solvated Parallel Displaced']

        # import pylab
        #
        # fig = pylab.figure()
        # figlegend = pylab.figure(figsize=(7.75, 0.6))
        # ax = fig.add_subplot(111)
        # for i in range(5):
        #     ax.plot(range(10), pylab.randn(10), color=colors[i], label=names[i])
        #
        # #lines = ax.plot(range(10), pylab.randn(10), range(10), pylab.randn(10),range(10), pylab.randn(10), range(10), pylab.randn(10), range(10), pylab.randn(10))
        # figlegend.legend(*ax.get_legend_handles_labels(), 'best', ncol=3)
        # fig.show()
        # figlegend.show()
        # fig.tight_layout()
        # figlegend.savefig('legend.pdf')
        #
        # plt.show()
        # exit()

        # names = ['Parallel Displaced (d=3.7)', 'Sandwiched (d=3.7)', 'Sandwiched (d=5.0)', 'Parallel Displaced (d=5.0)',
        #          'Solvated Parallel Displaced']
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
            print(i)
            fig = plt.figure(i)
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

            #plt.legend()#fontsize=12)

            plt.ylabel('Component Number Density \n (number/nm$^3$)', fontsize=18)
            plt.xlabel('Distance from pore center, r (nm)', fontsize=18)
            plt.axes().tick_params(labelsize=14)
            # plt.ylim([0, 0.6])
            plt.tight_layout()
            plt.savefig("%s_density.pdf" % regions[i])

        plt.show()

    if args.solvate:
        colors = ['xkcd:red', 'xkcd:green', 'blue', 'xkcd:yellow']
    else:
        colors = ['red', 'green', 'blue']

    if not args.load:

        R = Region(args.traj, args.gro, args.components, begin=args.begin, end=args.end)

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

        p_centers = physical.avg_pore_loc(npores, t.xyz[:, comp, :], box)

        for i, reg in enumerate(regions):

            print('Calculating number density of %s region' % reg)

            pos = p2p.restrict_atoms(t, reg)  # restrict trajectory to region

            p = duplicate(pos, t.unitcell_vectors)  # duplicate things periodically

            equil = 0
            density, r, bin_width = compdensity(pos, p_centers, equil, box, pores=npores, nbins=args.bins)

            results[i, :] = density

            plt.bar(r, density, bin_width, color=colors[i], alpha=0.6, label=reg)

        np.savez_compressed("regional_density", results=results, r=r, bw=bin_width, box=t.unitcell_vectors)
        print('Arrays saved as regional_density.npz')

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
