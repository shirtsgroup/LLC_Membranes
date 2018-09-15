#!/usr/bin/env python

import argparse
import numpy as np
from LLC_Membranes.analysis import disorder
import pickle
import tqdm
import os
import matplotlib.pyplot as plt
from scipy import stats


def initialize():

    parser = argparse.ArgumentParser(description='Crosslink LLC structure')  # allow input from user

    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Name of coordinate file')
    parser.add_argument('-t', '--traj', default='wiggle.trr', nargs='+', type=str, help='Name of trajectory to analyze')
    parser.add_argument('-b', '--begin', default=0, type=int, help='Start frame')
    parser.add_argument('-r', '--ref_atoms', nargs='+', default=['C', 'C1', 'C2', 'C3', 'C4', 'C5'], help='Name of '
                        'atoms that will be used to define head groups.')
    parser.add_argument('-pores', default=4, type=int, help='Number of pores in unit cell')
    parser.add_argument('-layers', default=20, type=int, help='Number of monomers per column')
    parser.add_argument('-out', default='disorder.png', type=str, help='')
    parser.add_argument('-fr', action="store_true", help='Force recompute. Even if trajectories pickle files exist,'
                                                         'redo the calculations')
    parser.add_argument('-pd', '--parallel_displaced', action="store_true", help='Specify if initial configuration is'
                                                                                 'parallel displaced')

    return parser


class Cdf(object):

    def __init__(self, data):
        """
        Generate an emperical cumulative distribution function from data
        :param data: x-values of data in no particular order
        """

        self.xs = np.array(sorted(data))
        N = float(len(self.xs))
        self.ys = np.arange(1, N + 1) / N

    def cdf(self, x):
        """
        Callable cumulative emperical distribution function
        :param x: array of x-values at which to evaluate cumulative emperical distribution function
        :return:
        """

        ndx = [np.argmin(np.abs(self.xs - x[i])) for i in range(x.size)]

        return self.ys[ndx]

    def random_sample(self, n=1):
        """
        :param n: number of random samples to draw (default=1)
        :return: random samples
        """

        return np.random.choice(self.xs, size=n, replace=True)


if __name__ == "__main__":

    args = initialize().parse_args()

    if os.path.isfile('%s/trajectories' % os.getcwd()) and not args.fr:
        print('Reloading Trajectories...', end='', flush=True)
        independent_trajectories = pickle.load(open("trajectories", "rb"))
        print('Done!')
    else:
        print('Analyzing Trajectories')
        independent_trajectories = []
        for i in tqdm.tqdm(range(len(args.traj)), unit='Trajectory'):
            independent_trajectories.append(disorder.System(args.traj[i], args.gro, args.ref_atoms, begin=args.begin,
                                                           pores=args.pores, layers=args.layers))
            independent_trajectories[i].z_deviation()
            independent_trajectories[i].r_deviation()
            independent_trajectories[i].theta_deviation(pd=args.parallel_displaced)

        with open('trajectories', 'wb') as f:
            pickle.dump(independent_trajectories, f)

    full_distribution = [[], [], []]  # z, r, theta
    sigmas = np.zeros([len(independent_trajectories), 3])
    for i in range(len(independent_trajectories)):
        full_distribution[0] += independent_trajectories[i].z_values.flatten().tolist()
        full_distribution[1] += independent_trajectories[i].r_values.flatten().tolist()
        full_distribution[2] += independent_trajectories[i].theta_values.flatten().tolist()
        # full_distribution[0] += independent_trajectories[i].z_values[-1, :].tolist()
        # full_distribution[1] += independent_trajectories[i].r_values[-1, :].tolist()
        # full_distribution[2] += independent_trajectories[i].theta_values[-1, ...].flatten().tolist()
        sigmas[i, :] = [independent_trajectories[i].dz, independent_trajectories[i].dr,
                        independent_trajectories[i].dtheta]

    full_distribution = [np.array(full_distribution[i]) for i in range(len(full_distribution))]
    full_distribution[2] -= full_distribution[2].mean()

    full_sigmas = [np.std(i) for i in full_distribution]  # standard deviation of full distribution

    means = np.zeros([len(independent_trajectories), 3])
    for i in range(means.shape[0]):
        means[i, 0] = independent_trajectories[i].z_values.mean()
        means[i, 1] = independent_trajectories[i].r_values.mean()
        means[i, 2] = independent_trajectories[i].theta_values.mean()

    means[:, 2] -= means[:, 2].mean()

    # For saving figures
    name = 'sandwiched'
    if args.parallel_displaced:
        name = 'offset'

    # Histograms of means of distributions
    fig, axmean = plt.subplots(1, 2, sharey=True)

    bins = 10
    # axmean[0].hist(means[:, 0], bins=bins)
    # axmean[0].plot([2.09e-17, 2.09e-17], [0, 13], '--', color='black', linewidth=2)  # parallel displaced
    # axmean[0].set_xlabel('$\mu_z$ (nm)', fontsize=16)
    # axmean[0].set_ylabel('Count', fontsize=16)
    # axmean[0].xaxis.set_tick_params(labelsize=12)
    # axmean[0].yaxis.set_tick_params(labelsize=12)

    axmean[0].hist(means[:, 1], bins=bins)
    axmean[0].plot([0.71, 0.71], [0, 13], '--', color='black', linewidth=2)  # parallel displaced
    axmean[0].set_xlabel('$\mu_r$ (nm)', fontsize=16)
    axmean[0].xaxis.set_tick_params(labelsize=12)
    axmean[0].yaxis.set_ticks(np.linspace(0, 8, 5))

    axmean[1].hist(means[:, 2], bins=bins)
    axmean[1].plot([-0.04, -0.04], [0, 12], '--', color='black', linewidth=2)  # parallel displaced
    axmean[1].set_xlabel('$\mu_\Theta$ (radians)', fontsize=16)
    axmean[1].xaxis.set_tick_params(labelsize=12)

    plt.ylim(0, 9)
    plt.tick_params(labelsize=14)
    plt.tight_layout()
    plt.savefig(
        '/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/%s_ensemble_means.png' % name)

    # Histograms of standard deviations of distributions
    fig, ax0 = plt.subplots(1, 3, sharey=True)

    print('Standard Deviation of standard deviations:')
    print('Sigma_z = %.3f' % np.std(sigmas[:, 0]))
    print('Sigma_r = %.3f' % np.std(sigmas[:, 1]))
    print('Sigma_theta = %.3f' % np.std(sigmas[:, 2]))

    bins = 10
    ax0[0].hist(sigmas[:, 0], bins=bins)
    ax0[0].plot([0.1628, 0.1628], [0, 12], '--', color='black', linewidth=2)
    ax0[0].set_xlabel('$\sigma_z$ (nm)', fontsize=16)
    ax0[0].set_ylabel('Count', fontsize=16)
    ax0[0].xaxis.set_tick_params(labelsize=12)
    ax0[0].yaxis.set_tick_params(labelsize=12)

    ax0[1].hist(sigmas[:, 1], bins=bins)
    ax0[1].plot([0.2060, 0.2060], [0, 12], '--', color='black', linewidth=2)
    ax0[1].set_xlabel('$\sigma_r$ (nm)', fontsize=16)
    ax0[1].xaxis.set_tick_params(labelsize=12)

    ax0[2].hist(sigmas[:, 2], bins=bins)
    ax0[2].plot([0.2429, 0.2429], [0, 12], '--', color='black', linewidth=2)
    ax0[2].set_xlabel('$\sigma_\Theta$ (radians)', fontsize=16)
    ax0[2].xaxis.set_tick_params(labelsize=12)

    plt.ylim(0, 11.5)
    plt.tick_params(labelsize=14)
    plt.tight_layout()

    plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/%s_ensemble_stds.png' % name)

    # TEST INDIVIDUAL CENTER OF MASSES
    z_cdf = Cdf(full_distribution[0])
    r_cdf = Cdf(full_distribution[1])
    theta_cdf = Cdf(full_distribution[2])

    nframes = 51
    z = np.zeros([len(independent_trajectories), 400])
    r = np.zeros_like(z)
    theta = np.zeros_like(z)

    for i in range(len(independent_trajectories)):
        z[i, :] = independent_trajectories[i].z_values[-1, :]  # value at last frame for each center of mass
        r[i, :] = independent_trajectories[i].r_values[-1, :]  # value at last frame for each center of mass
        theta[i, :] = independent_trajectories[i].theta_values[-1, ...].flatten()  # value at last frame for each center of mass

    N = len(independent_trajectories)
    trials = 1000
    boot = np.zeros([trials, N, 3])
    x_ecdf = np.zeros([N, 3])
    for i in range(N):
        x_ecdf[i, 0] = z_cdf.xs[np.argmin(np.abs(z_cdf.ys - (i + 1) / N))]
        x_ecdf[i, 1] = r_cdf.xs[np.argmin(np.abs(r_cdf.ys - (i + 1) / N))]
        x_ecdf[i, 2] = theta_cdf.xs[np.argmin(np.abs(theta_cdf.ys - (i + 1) / N))]

    all_ecdfs = np.zeros([400, N, 3])
    for i in range(all_ecdfs.shape[0]):
        all_ecdfs[i, :, 0] = np.array(sorted(z[:, i]))
        all_ecdfs[i, :, 1] = np.array(sorted(r[:, i]))
        all_ecdfs[i, :, 2] = np.array(sorted(theta[:, i]))

    # plt.hist(all_ecdfs.mean(axis=1)[:, 1], bins=50)
    # plt.show()
    # ys = np.arange(1, N + 1) / N
    # choices = np.random.choice(all_ecdfs.shape[0], size=10)
    # print(choices)
    # for i in range(len(choices)):
    #     plt.plot(all_ecdfs[choices[i], :, 0], ys)
    # plt.show()
    # exit()

    # print([np.mean(full_distribution[i]) for i in range(3)])
    # print(mean_ecdf.mean(axis=0))
    #
    # print([np.std(full_distribution[i]) for i in range(3)])
    # print(mean_ecdf.std(axis=0))
    # exit()

    confidence = 95  # percent confidence interval
    lower_confidence = (100 - confidence) / 2
    upper_confidence = 100 - lower_confidence

    mean_ecdf = all_ecdfs.mean(axis=0)

    error_ecdf = np.zeros([2, N, 3])
    error_ecdf[0, :, 0] = np.abs(np.percentile(all_ecdfs[..., 0], lower_confidence, axis=0) - mean_ecdf[:, 0])  # 2.5 percent of data below this value
    error_ecdf[1, :, 0] = np.percentile(all_ecdfs[..., 0], upper_confidence, axis=0) - mean_ecdf[:, 0]  # 97.5 percent of data below this value
    error_ecdf[0, :, 1] = np.abs(np.percentile(all_ecdfs[..., 1], lower_confidence, axis=0) - mean_ecdf[:, 1])
    error_ecdf[1, :, 1] = np.percentile(all_ecdfs[..., 1], upper_confidence, axis=0) - mean_ecdf[:, 1]
    error_ecdf[0, :, 2] = np.abs(np.percentile(all_ecdfs[..., 2], lower_confidence, axis=0) - mean_ecdf[:, 2])
    error_ecdf[1, :, 2] = np.percentile(all_ecdfs[..., 2], upper_confidence, axis=0) - mean_ecdf[:, 2]

    ys = np.arange(1, N + 1) / N
    shift = 0.0025
    #
    # plt.errorbar(mean_ecdf[:, 0], ys, xerr=error_ecdf[..., 0])
    # plt.plot(z_cdf.xs, z_cdf.ys, color='black')
    # plt.show()
    # exit()

    # find where y = [1/40, 2/40, 3/40 ...]
    for i in range(trials):
        boot[i, :, 0] = np.array(sorted(z_cdf.random_sample(N)))
        boot[i, :, 1] = np.array(sorted(r_cdf.random_sample(N)))
        boot[i, :, 2] = np.array(sorted(theta_cdf.random_sample(N)))

    error = np.zeros([2, N, 3])
    error[0, :, 0] = np.abs(np.percentile(boot[..., 0], lower_confidence, axis=0) - x_ecdf[:, 0])  # 2.5 percent of data below this value
    error[1, :, 0] = np.percentile(boot[..., 0], upper_confidence, axis=0) - x_ecdf[:, 0]  # 97.5 percent of data below this value
    error[0, :, 1] = np.abs(np.percentile(boot[..., 1], lower_confidence, axis=0) - x_ecdf[:, 1])
    error[1, :, 1] = np.percentile(boot[..., 1], upper_confidence, axis=0) - x_ecdf[:, 1]
    error[0, :, 2] = np.abs(np.percentile(boot[..., 2], lower_confidence, axis=0) - x_ecdf[:, 2])
    error[1, :, 2] = np.percentile(boot[..., 2], upper_confidence, axis=0) - x_ecdf[:, 2]

    fig, axboot = plt.subplots(1, 3, figsize=(11, 4), sharey=False)
    axboot[0].errorbar(x_ecdf[:, 0], z_cdf.cdf(x_ecdf[:, 0]) - shift, xerr=error[..., 0], ecolor='black', elinewidth=1,
                       color='black', linewidth=0)
    axboot[0].errorbar(mean_ecdf[:, 0], ys + shift, xerr=error_ecdf[..., 0], ecolor='red', linewidth=2, elinewidth=1, color='red')
    axboot[0].plot(z_cdf.xs, z_cdf.ys, linewidth=2, color='black')
    axboot[0].set_ylabel('Cumulative Probability', fontsize=14)
    axboot[0].set_xlim(-0.65, 0.65)
    axboot[0].set_xlabel('z-deviation (nm)', fontsize=14)
    axboot[0].xaxis.set_tick_params(labelsize=12)
    axboot[0].yaxis.set_tick_params(labelsize=12)

    axboot[1].errorbar(x_ecdf[:, 1], r_cdf.cdf(x_ecdf[:, 1]) - shift, xerr=error[..., 1], ecolor='black', elinewidth=1,
                       color='black', linewidth=0)
    axboot[1].errorbar(mean_ecdf[:, 1], ys + shift, xerr=error_ecdf[..., 1], ecolor='red', linewidth=2, elinewidth=1, color='red')
    axboot[1].plot(r_cdf.xs, r_cdf.ys, linewidth=2, color='black')
    axboot[1].set_xlim(0, 1.25)
    axboot[1].set_xlabel('r-deviation (nm)', fontsize=14)
    axboot[1].xaxis.set_tick_params(labelsize=12)
    axboot[1].yaxis.set_tick_params(labelsize=12)

    axboot[2].errorbar(x_ecdf[:, 2], theta_cdf.cdf(x_ecdf[:, 2]) - shift, xerr=error[..., 2], ecolor='black', elinewidth=1,
                       color='black', linewidth=0)
    axboot[2].errorbar(mean_ecdf[:, 2], ys + shift, xerr=error_ecdf[..., 2], ecolor='red', linewidth=2, elinewidth=1, color='red')
    axboot[2].plot(theta_cdf.xs, theta_cdf.ys, linewidth=2, color='black')
    axboot[2].set_xlim(-np.pi/2, np.pi/2)
    axboot[2].set_xlabel('$\Theta$-deviation (radians)', fontsize=14)
    axboot[2].xaxis.set_tick_params(labelsize=12)
    axboot[2].yaxis.set_tick_params(labelsize=12)
    axboot[2].set_xticklabels(['-$\pi$/2', '-$\pi$/4', '0', '$\pi$/4', '$\pi$/2'])

    nplot = 3  # number of ecdf's to plot
    random_ecdfs = np.random.choice(400, size=(nplot, 3))
    print(random_ecdfs)

    if name == 'sandwiched':
        z_ecdfs = [394, 315, 211, 218, 240]
        r_ecdfs = [79, 53, 37, 160, 309]
        theta_ecdfs = [140, 39, 107, 157, 287]
    else:
        # parallel displaced
        z_ecdfs = [204, 270, 92, 109, 101]
        r_ecdfs = [244, 70, 294, 55, 22]  # 22, 189
        theta_ecdfs = [140, 39, 107, 157, 287]

    # for i in z_ecdfs:
    for i in random_ecdfs[:, 0]:
    # i = np.random.choice(400)
    # print(i)
        x_emp = np.array(sorted(z[:, i]))
        n = float(len(x_emp))
        ys = np.arange(1, n + 1) / n
        axboot[0].plot(x_emp, ys, linewidth=2)

    #for i in r_ecdfs:
    for i in random_ecdfs[:, 1]:
    # i = np.random.choice(400)
    # print(i)
        x_emp = np.array(sorted(r[:, i]))
        n = float(len(x_emp))
        ys = np.arange(1, n + 1) / n
        axboot[1].plot(x_emp, ys, linewidth=2)

    #for i in theta_ecdfs:
    for i in random_ecdfs[:, 2]:
        x_emp = np.array(sorted(theta[:, i]))
        n = float(len(x_emp))
        ys = np.arange(1, n + 1) / n
        axboot[2].plot(x_emp, ys, linewidth=2)

    plt.tight_layout()
    plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/%s_ecdfs.png' % name)
    # plt.show()
    # exit()

    # pz = []
    # pr = []
    # ptheta = []
    # for i in tqdm.tqdm(range(independent_trajectories[0].z_values.shape[1])):
    #     pz.append(stats.kstest(z[:, i], lambda x: z_cdf.cdf(x))[1])
    #     # if 0.45 < pz[i] < 0.55:
    #     #     x_emp = np.array(sorted(z[:, i]))
    #     #     N = float(len(x_emp))
    #     #     ys = np.arange(1, N + 1) / N
    #     #     plt.plot(x_emp, ys)
    #     #     plt.plot(x_emp, z_cdf.cdf(x_emp))
    #     #     plt.xlabel('x')
    #     #     plt.ylabel('cumulated probability')
    #     #     plt.show()
    #     #     exit()
    #     pr.append(stats.kstest(r[:, i], lambda x: r_cdf.cdf(x))[1])
    #     ptheta.append(stats.kstest(theta[:, i], lambda x: theta_cdf.cdf(x))[1])

    # plt.figure()
    # for i in range(10):
    #     outline = np.zeros([z.shape[0] * 2 + 2, 2])
    #     plt.hist(z[:, i])
    # plt.xlabel('z-deviation (nm)')
    # plt.ylabel('Frequency')
    #
    # plt.show()

    # fig, ax0 = plt.subplots(1, 3, figsize=(12, 5))
    # ax0[0].hist(pz, bins=20)
    # ax0[0].set_xlabel('p-value')
    # ax0[0].set_ylabel('Count')
    # ax0[0].set_title('z')
    # ax0[1].hist(pr, bins=20)
    # ax0[1].set_xlabel('p-value')
    # ax0[1].set_title('r')
    # ax0[2].hist(ptheta, bins=20)
    # ax0[2].set_xlabel('p-value')
    # ax0[2].set_title('$\Theta$')
    # plt.show()

    #
    # for i in range(400):
    #     for j in range(len(independent_trajectories))
    #     independent_trajectories[]

    # print(stats.kstest(independent_trajectories[0].z_values.flatten(), lambda x: z_cdf.cdf(x)))
    # plt.plot(z_cdf.xs, z_cdf.ys)
    # x_emp = sorted(independent_trajectories[0].z_values.flatten().tolist())
    # N = float(len(x_emp))
    # ys = np.arange(1, N + 1) / N
    # plt.plot(x_emp, ys)
    # plt.show()
    # exit()

    # compare each frame of each individual distribution to full distribution
    # p = [[], [], []]
    # for i in tqdm.tqdm(range(len(independent_trajectories))):
    #     ecdf = Cdf(independent_trajectories[i].z_values.flatten())
    #     plt.plot(ecdf.xs, ecdf.ys)
        # for t in range(independent_trajectories[i].t.n_frames):
        #     p[0].append(stats.kstest(independent_trajectories[0].z_values[t, :], lambda x: z_cdf.cdf(x))[1])
            #p[0].append(stats.ks_2samp(independent_trajectories[i].z_values[t, :], full_distribution[0])[1])
            #p[1].append(stats.ks_2samp(independent_trajectories[i].r_values[t, :], full_distribution[1]))
            #p[2].append(stats.ks_2samp(independent_trajectories[i].z_values[t, :], full_distribution[0])[1])

    # plt.plot(z_cdf.xs, z_cdf.ys, linewidth=2, label='All data combined', color='black')
    # plt.legend()
    # plt.xlabel('x')
    # plt.ylabel('Cumulated distribution')
    # plt.show()
    # exit()
    # plt.hist(p[0], bins=50)
    # plt.show()
    # exit()
    # compare distribution of all data combined from each system to full distribution. Note: array might need flattening
    # p = [[], [], []]
    # for i in range(len(independent_trajectories)):
    #     p[0].append(stats.ks_2samp(independent_trajectories[i].z_values.tolist(), full_distribution[0])[1])
    #     p[1].append(stats.ks_2samp(independent_trajectories[i].r_values.tolist(), full_distribution[1])[1])
    #     p[2].append(stats.ks_2samp(independent_trajectories[i].theta_values.tolist(), full_distribution[2])[1])

    # Compare individual distributions against each other individual distribution
    # p2d = np.zeros([len(independent_trajectories), len(independent_trajectories), 3])
    # for i in range(len(independent_trajectories)):
    #     for j in range(len(independent_trajectories)):
    #         if i != j:
    #             p2d[i, j, 0] = stats.ks_2samp(independent_trajectories[i].z_values.tolist(),
    #                                           independent_trajectories[j].z_values.tolist())[1]
    #             p2d[i, j, 1] = stats.ks_2samp(independent_trajectories[i].r_values.tolist(),
    #                                           independent_trajectories[j].r_values.tolist())[1]
    #             p2d[i, j, 2] = stats.ks_2samp(independent_trajectories[i].theta_values.tolist(),
    #                                           independent_trajectories[j].theta_values.tolist())[1]

    # PLOT FULL DISTRIBUTIONS
    fig, ax = plt.subplots(1, 3, figsize=(11, 4), sharey=False)

    bins = 50
    ax[0].hist(full_distribution[0], bins=bins, normed=True, range=(-0.6, 0.6))
    #ax[0].set_title('$\sigma=%.4f nm$' % np.std(full_distribution[0]))
    ax[0].set_xlabel('z-deviation (nm)', fontsize=14)
    ax[0].set_ylabel('Frequency', fontsize=14)
    ax[0].xaxis.set_tick_params(labelsize=12)
    ax[0].yaxis.set_tick_params(labelsize=12)

    ax[1].hist(full_distribution[1], bins=bins, normed=True)
    #ax[1].set_title('$\sigma=%.4f nm$' % np.std(full_distribution[1]))
    ax[1].set_xlabel('r-deviation (nm)', fontsize=14)
    ax[1].xaxis.set_tick_params(labelsize=12)
    ax[1].yaxis.set_tick_params(labelsize=12)
    ax[1].yaxis.set_ticks(np.linspace(0, 2, 5))

    ax[2].hist(full_distribution[2], bins=bins, normed=True, range=(-np.pi/2, np.pi/2))
    #ax[2].set_title('$\sigma=%.4f nm$' % np.std(full_distribution[2]))
    ax[2].set_xlabel('$\Theta$-deviation (radians)', fontsize=14)
    ax[2].xaxis.set_tick_params(labelsize=12)
    ax[2].yaxis.set_tick_params(labelsize=12)
    ax[2].set_xticklabels(['-$\pi$/2', '-$\pi$/4', '0', '$\pi$/4', '$\pi$/2'])
    ax[2].yaxis.set_ticks(np.linspace(0, 1.5, 4))

    plt.tight_layout()
    plt.savefig('/home/bcoscia/PycharmProjects/LLC_Membranes/Ben_Manuscripts/structure_paper/figures/%s_ensemble_pooled.png' % name)
    plt.show()
    exit()
    fig, ax2 = plt.subplots(1, 3)
    ax2[0].hist(p[0], bins=bins)
    ax2[0].set_title('z-direction')
    ax2[0].set_xlabel('p-value')
    ax2[0].set_ylabel('Frequency')

    ax2[1].hist(p[1], bins=bins)
    ax2[1].set_title('r-direction')
    ax2[1].set_xlabel('p-value')

    ax2[2].hist(p[2], bins=bins)
    ax2[2].set_title('$\Theta$-direction')
    ax2[2].set_xlabel('p-value')

    # fig, ax3 = plt.subplots(1, 3)
    # im1 = ax3[0].imshow(p2d[..., 0])
    # ax3[0].set_title('z-direction')
    #
    # im2 = ax3[1].imshow(p2d[..., 1])
    # ax3[1].set_title('r-direction')
    #
    # im3 = ax3[2].imshow(p2d[..., 2])
    # ax3[2].set_title('$\Theta$-direction')
    #
    # fig.colorbar(im1, ax=ax3[0])
    # fig.colorbar(im2, ax=ax3[1])
    # fig.colorbar(im3, ax=ax3[2])

    plt.tight_layout()
    plt.show()