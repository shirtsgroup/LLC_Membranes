#!/usr/bin/env python

# Good explanation/example of kolmogorov-smirnov test : https://onlinecourses.science.psu.edu/stat414/node/323/

import argparse
import numpy as np
import mdtraj as md
from LLC_Membranes.analysis import residues, p2p
from LLC_Membranes.timeseries import correlation
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from statsmodels.tsa.stattools import adfuller


def initialize():

    parser = argparse.ArgumentParser(description='Test if particles undergo brownian motion')

    parser.add_argument('-t', '--traj', default='npt_nojump.xtc', type=str, help='Name of trajectory file (.xtc or '
                        '.trr). NOTE: Use unwrapped coordinates')
    parser.add_argument('-g', '--gro', default='initial.gro', type=str, help='Coordinate file')
    parser.add_argument('-r', '--residue', default='HOH', type=str, help='Name of residue whose center of mass you want'
                                                                     'to track')
    parser.add_argument('-l', '--lag', default=1, type=int, help='Number of increments between calculation of deviation'
                                                                'from previous position')
    parser.add_argument('--restrict_to_pores', action="store_true", help='Only analyze water molecules that are '
                                                                         'within HII phase pores')
    parser.add_argument('-radius', '--radius', default=1, type=float, help='Pore radius (nm) when restricting COMs to pore'
                                                                      ' centers')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Frame at which to start calculations')
    parser.add_argument('-e', '--end', default=-1, type=int, help='Last frame to include in calculations')

    return parser


def center_of_mass(pos, mass):
    """
    Calculate center of mass of groups of atoms
    :param pos: atomic positions of resiude whose center of mass is to be calculated
    :param mass: mass of each atom in the residue
    :return: trajectory of center of mass coordinates
    """

    n = len(mass)  # number of atoms per residue
    nres = pos.shape[1] // n  # number of residues
    com = np.zeros([pos.shape[0], nres, 3])  # center of mass array

    for f in range(pos.shape[0]):
        for g in range(nres):  # calculate center of mass of each residue
                start = g * n
                end = (g + 1) * n
                w = (pos[f, start:end, :].T * mass).T
                com[f, g, :] = np.sum(w, axis=0) / sum(mass)

    return com


def gaussian(x, sigma, mean):
    """

    :param x: independent variable, (float or np.array)
    :param sigma: standard deviation
    :param mean: mean
    :return value(s) of gaussian at x value(s)
    """

    return (1 / np.sqrt(2 * np.pi * sigma ** 2)) * np.exp(-(x - mean)**2 / (2 * sigma ** 2))


def restrict_to_pore(t, com, r, npores=4, ref_atoms=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], ref_frame=-1):
    """

    :param t: mdtraj trajectory object for full system
    :param com: center of mass trajectory of residue of interest
    :param r: distance from pore center to be considered inside the pore region
    :param npores: number of pores in HII system
    :param ref_atoms: reference atoms used to determine pore centers
    :param ref_frame: reference frame for determining pore centers. This is an approximation since the system
    changes dynamically
    :return: center of mass trajectory excluding COM's outside of the pore
    """

    ref = [a.index for a in t.topology.atoms if a.name in ref_atoms]
    pcenters = p2p.avg_pore_loc(npores, t.xyz[:, ref, :], 0)

    cut = []
    for p in range(npores):
        d = np.linalg.norm(com[ref_frame, :, :2] - pcenters[:, p, ref_frame], axis=1)
        cut += np.where(d < r)[0].tolist()

    #opposite = [i for i in range(com.shape[1]) if i not in cut]
    #return com[:, opposite, :]

    return com[:, cut, :]


if __name__ == "__main__":

    args = initialize().parse_args()

    # load up trajectory
    t = md.load(args.traj, top=args.gro)[args.begin:args.end]  # load trajectory
    nframes = t.n_frames
    dt = t.time[-1] / (nframes - 1)

    # restrict trajectory to atoms that make up residue of interest
    res = residues.Residue(args.residue)  # residue object with all of its relevant properties
    if res.resname == 'SOL':  # mdtraj workaround
        res.resname = 'HOH'

    keep = [a.index for a in t.topology.atoms if a.residue.name == res.resname]  # keep only atoms that are part of res

    nres = len(keep) / res.natoms

    pos = t.xyz[:, keep, :]

    com = center_of_mass(pos, res.masses)

    if args.restrict_to_pores:
        com = restrict_to_pore(t, com, args.radius)

    # -- plot random molecule COM traces -- #
    # n_traces = 5
    # traj = np.random.choice(com.shape[1], size=n_traces, replace=False)
    # for i in range(n_traces):
    #     plt.plot(t.time/1000, com[:, traj[i], 2], linewidth=2)
    # plt.xlabel('Time (ns)', fontsize=14)
    # plt.ylabel('$z$-coordinate (nm)', fontsize=14)
    # plt.gcf().get_axes()[0].tick_params(labelsize=14)
    # plt.tight_layout()
    # plt.show()
    # exit()

    # calculate deviation in x, y and z at each frame including all residue molecules at each frame
    # in other words, cacluate the args.lag-th order difference equation
    deviation = np.zeros([com.shape[0] - args.lag, com.shape[1], 3])
    for i in range(args.lag, nframes):
            deviation[i - args.lag, ...] = com[i, ...] - com[i - args.lag, ...]

    # -- Augmented Dickey-Fuller for each individual molecule -- #
    # adf = []
    # for i in range(com.shape[1]):
    #     p = adfuller(com[:, i, 2])[1]
    #     if p < 0.05:
    #         plt.plot(t.time/1000, com[:, i, 2])
    #         plt.show()
    #     adf.append(p)
    #
    # plt.hist(adf, bins=20)
    # plt.ylabel('Count')
    # plt.xlabel('$p$-value')
    # plt.tight_layout()
    # plt.savefig('water_p_hist.pdf')
    # plt.show()
    # exit()

    # -- Augmented Dickey-Fuller test using all data -- #
    # result = adfuller(com[..., 2].flatten())  # want p-value as low as possible if stationary
    # print('ADF Statistic: %f' % result[0])
    # print('p-value: %f' % result[1])
    # print('Critical Values:')
    # for key, value in result[4].items():
    #     print('\t%s: %.3f' % (key, value))

    # -- Autocorrelation functions for each molecule -- #
    for i in range(com.shape[1]):
        print(i)
        fig, ax = plt.subplots(1, 2, figsize=(8, 4))
        ax[0].plot(t.time/1000, com[:, i, 2], linewidth=2)
        ax[0].set_xlabel('Time (ns)', fontsize=14)
        ax[0].set_ylabel('$z$-coordinate (nm)', fontsize=14)
        ax[0].set_title('Time Series', fontsize=14)
        ax[0].xaxis.set_tick_params(labelsize=12)
        ax[0].yaxis.set_tick_params(labelsize=12)

        # ax[1].plot(t.time, correlation.acf(com[:, i, 2]))
        ax[1].plot(t.time[:-(args.lag + 1)]/1000, correlation.acf(deviation[:, i, 2]), linewidth=2)
        ax[1].set_xlabel('Lag (ns)', fontsize=14)
        ax[1].set_title('Autocorrelation Function', fontsize=14)
        plt.tight_layout()
        if i == 6:
            plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/transport/figures/eth_autocorrelation.pdf')
        # ax[1].set_xlim(0, 200)
        plt.show()
    exit()

    # # -- Autocovariance using all increments -- #
    # fig, ax = plt.subplots(1, 3, figsize=(12, 5))
    # ax[0].plot(t.time/1000, com[:, 1, 2])
    # ax[0].set_xlabel('time (ns)')
    # ax[0].set_title('$z$-coordinate vs. time')
    # ax[1].plot(t.time[:-args.lag]/1000, correlation.autocov(deviation[..., 2].T))
    # ax[1].set_xlabel('time (ns)')
    # ax[1].set_title('$\gamma_j$')
    # ax[2].plot(t.time[:-args.lag]/1000, np.mean(deviation[..., 2], axis=1))
    # ax[2].set_title('$\mu$')
    # ax[2].set_xlabel('time (ns)')
    # plt.tight_layout()
    # plt.show()
    # exit()

    # -- Autocorrelation functions for a fixed lag. Stats calculated by bootstrapping -- #
    # fig, ax = plt.subplots()
    #
    # nframes = int(com.shape[0] / args.lag)
    # if args.lag > 1:
    #     nframes += 1
    # nmolecules = com.shape[1]
    # acfs = np.zeros([nframes, nmolecules])  # autocorrelation function for each molecule
    # for i in range(nmolecules):
    #     acfs[:, i] = correlation.acf(com[::args.lag, i, 2])

    # plot a bunch of ACFs
    # show = int(0.5 * nframes)
    # for i in range(50):
    #     plt.plot(t.time[:int(0.5*t.time.size):args.lag], acfs[:show, i])

    # nboot = 200
    # acf_boot = np.zeros([nboot, nframes])
    #
    # # for i in range(nframes):
    # #     acf_boot[:, i] = acfs[i, np.random.choice(nmolecules, size=nboot, replace=True)]
    #
    # for i in range(nboot):
    #     choices = np.random.choice(nmolecules, size=nmolecules, replace=True)  # randomly choose water trajectories with replacement
    #     acf_boot[i, :] = acfs[:, choices].mean(axis=1)
    #
    # error_boot = np.zeros([2, nframes])
    # error_boot[0, :] = np.abs(np.percentile(acf_boot, 2.5, axis=0) - np.mean(acf_boot, axis=0)) # 2.5 percent of data below this value
    # error_boot[1, :] = np.percentile(acf_boot, 97.5, axis=0) - np.mean(acf_boot, axis=0) # 97.5 percent of data below this value
    #
    # bootstrapped = np.mean(acf_boot, axis=0)
    #
    # show = int(0.5 * t.time.size)
    # ax.errorbar(t.time[:show] / 1000, bootstrapped[:show], yerr=error_boot[:, :show], ecolor='black', elinewidth=1,
    #             color='blue', linewidth=2, errorevery=5)
    # ax.set_xlabel('time (ns)')
    # plt.tight_layout()
    # plt.savefig('acf_eth_HII.pdf')
    # plt.show()

    # -- Plot distribution of frame-to-frame deviations -- #
    fig, ax = plt.subplots(1, 3, sharey=True, figsize=(12, 5))
    direction = {0: 'x', 1: 'y', 2: 'z'}
    sigmas = [np.std(deviation[..., i].flatten()) for i in range(3)]

    for i in range(3):
        n, bins, _ = ax[i].hist(deviation[..., i].flatten(), bins=50, normed=True)
        ax[i].set_title('$\sigma_%s$ = %.3f' % (direction[i], sigmas[i]))
        # fit gaussian to histogram
        p = [sigmas[i], 0]  # initial guess at sigma and mean
        x = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]  # centers of bins
        popt, pcov = curve_fit(gaussian, x, n, p0=p)
        ax[i].plot(x, gaussian(x, popt[0], popt[1]), linewidth=2)
        ax[i].set_xlabel('Deviation (nm)', fontsize=14)
        plt.gcf().get_axes()[i].tick_params(labelsize=14)

    ax[0].set_ylabel('Probability', fontsize=14)


    plt.show()
    exit()

    # test null hypothesis (that incremental steps are normally distributed) for individual water trajectories
    # plt.figure()
    # plt.title('Distribution of $p$-values from individual water trajectories')
    fig, ax = plt.subplots(1, 3, sharey=True)

    means = []
    for d in range(3):
        p = []

        for i in range(deviation.shape[1]):
            pvalue = stats.kstest(deviation[:, i, d], 'norm', args=(0, np.std(deviation[:, i, d])))[1]
            p.append(pvalue)
            means.append(np.mean(deviation[:, i, d]))

        if d == 0:
            plt.figure()
            plt.hist(p, bins=50)
            plt.xlabel('$p$-value', fontsize=16)
            plt.ylabel('Number of trajectories', fontsize=16)
            plt.tight_layout()
            plt.savefig('p-distribution_lag_%s' % args.lag)

        ax[d].set_title('%s-direction' % direction[d])
        ax[d].hist(p, bins=50)
        ax[d].set_xlabel('$p$-value', fontsize=16)
        plt.tight_layout()

    ax[0].set_ylabel('Number of trajectories', fontsize=16)

    # hand calculated kolmogorov-smirnov test
    plt.figure()
    xs = sorted(deviation.flatten())  # sort all deviation values from lowest to highest
    N = float(len(xs))  # number of values in xs
    ys = np.arange(1, N + 1) / N  # empirical cumulative distribution function
    cdf = stats.norm(scale=np.std(xs)).cdf(xs)  # normal cumulative distribution function at each point
    plt.plot(xs, ys, label='ecdf')
    plt.plot(xs, cdf, label='cdf')
    plt.legend()

    diff1 = np.abs(ys - cdf)
    diff2 = np.abs(np.concatenate(([0], ys))[:-1] - cdf)

    D = max(np.amax(diff1), np.amax(diff2))

    print('D statistic: %.4f' % D)
    print('Critical Value: %.4f' % (
                1.35810 / np.sqrt(N)))  # http://www.real-statistics.com/statistics-tables/kolmogorov-smirnov-table/
    print('P-value : %.4f' % (stats.distributions.kstwobign.sf(
        D * np.sqrt(N))))  # https://github.com/scipy/scipy/blob/v0.14.0/scipy/stats/stats.py#L3307

    # ks test for all data
    data = deviation.flatten()
    print(stats.kstest(data, 'norm', args=(0, np.std(data))))
    plt.show()
