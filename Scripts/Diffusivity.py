#!/usr/bin/bash

import argparse
import numpy as np
from Get_Positions import get_positions as gp
import matplotlib.pyplot as plt
import Poly_fit


def initialize():

    parser = argparse.ArgumentParser(description = 'Run Cylindricity script')
    parser.add_argument('-i', '--input', default='wiggle_traj.gro', help = 'Path to input file')
    parser.add_argument('-c', '--comp', default='NA', help='Name of ion(s) being used to calculate ionic conductivity')
    parser.add_argument('-b', '--nboot', default=200, help='Number of bootstrap trials to be run')
    parser.add_argument('-f', '--frontfrac', default=0.05, help='Where to start fitting line for diffusivity calc')
    parser.add_argument('-F', '--fracshow', default=.2, help='Percent of graph to show, also where to stop fitting line'
                                                             'during diffusivity calculation')
    parser.add_argument('-d', '--dim', default=1, help='Number of dimensions in which diffusion is being calculated')
    parser.add_argument('-l', '--LC_type', default='HII', help='Type of liquid crystal being studied')
    parser.add_argument('-s', '--solv', default='no', help='Is the system solvated or not?')
    parser.add_argument('-m', '--nMC', default=1000, help='Number of Monte Carlo trials to estimate error in D')
    args = parser.parse_args()

    return args


def autocorrFFT(x):
    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   # now we have the autocorrelation in convention B
    n= N*np.ones(N)-np.arange(0, N)  # divide res(m) by (N-m)
    return res/n  # this is the autocorrelation in convention A


def msd_fft(r):
    N = len(r)
    D = np.square(r).sum(axis=1)
    D = np.append(D, 0)
    S2 = sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q = 2*D.sum()
    S1 = np.zeros(N)
    for m in range(N):
      Q = Q-D[m-1]-D[N-m]
      S1[m] = Q/(N-m)
    return S1-2*S2


def msd(x, nT, no_comp):
    MSD = np.zeros([nT], dtype=float)
    MSDs = np.zeros([nT, no_comp], dtype=float)  # a set of MSDs per particle
    for n in range(0, no_comp):
        MSDs[:, n] = msd_fft(x[:, n, :].T)
        MSD += MSDs[:, n]
    MSD /= no_comp
    return MSD, MSDs


def bootstrap(nT, Nbootstraps, no_comp, MSDs, MSD):
    eMSDs = np.zeros([nT, Nbootstraps], dtype=float)  # a set of MSDs per particle (time step?)
    # Now, bootstrap over number of particles, assuming that the particles are sufficiently independent
    for b in range(0, Nbootstraps):
        print 'bootstrap trial number %s' % (b + 1)
        indices = np.random.randint(0, no_comp, no_comp)  # randomly select N of the particles with replacement
        # ^ makes a list of length Nparticles each with a random number from 0 to  Nparticles
        for n in range(0, no_comp):
            eMSDs[:, b] += MSDs[:, indices[n]]  # for this bootstrap trial b, add the MSDs of a randomly selected particle
            # to the particle at each time step
        eMSDs[:, b] /= no_comp  # Divide every timestep by Nparticles -- average the MSDs

    limits = np.zeros([2, nT], dtype=float)  # a set of MSDs per particle
    # now, let's determine a 95\% error bound for each tau (out of the
    # Nbootstrapped MSD's, use that for the error bars
    for t in range(0, nT):
        limits[0, t] = np.abs(np.percentile(eMSDs[t, :], 2.5) - MSD[t])
        limits[1, t] = np.abs(np.percentile(eMSDs[t, :], 97.5) - MSD[t])

    return limits


def d_error(nMC, startfit, endMSD, nT, limits, times, MSD, d):

    # Monte Carlo Method
    # nMC = 1000  # number of Monte Carlo trajectories
    tot_pts = (int(endMSD) - int(startfit))
    # pts = np.zeros([nMC, tot_pts])
    # for i in range(nMC):
    #     pts[i, :] = limits[0, startfit:endMSD]*np.random.randn(1, tot_pts) + MSD[startfit:endMSD]
    #
    # slopes = np.zeros([nMC])
    # for i in range(nMC):
    #     slopes[i] = Poly_fit.poly_fit(times[startfit:endMSD], pts[i, :], 1)[4][1]
    #
    # avg_D = np.mean(slopes)/(2*d*1000000)  # take average and convert to m^2/s
    # std_D = np.std(slopes)/(2*d*1000000)
    # print 'Error in slope calculated'

    # Weight least squares -- first generate a weight matrix
    W = np.zeros((tot_pts, tot_pts))
    for i in range(tot_pts):
        W[i, i] = 1/((limits[0, i + startfit])**2)

    y_fit, _, slope_error, _, A = Poly_fit.poly_fit(times[startfit:endMSD], MSD[startfit:endMSD], 1, W)

    return A[1]/(2*d*1000000), slope_error/(2*d*1000000)


def dconst(x, nT, Nbootstraps, frontfrac, fracshow, d, dt, nMC):

    no_comp = len(x[0, :, 0])
    print x.shape[0]
    print x.shape[1]
    print x.shape[2]
    exit()
    MSD, MSDs = msd(x, nT, no_comp)
    print 'MSDs done'
    limits = bootstrap(nT, Nbootstraps, no_comp, MSDs, MSD)
    print 'Bootstrapping Done'
    endMSD = int(np.floor(nT*fracshow))
    startMSD = int(np.floor(nT*frontfrac))
    # coeff = Poly_fit.poly_fit(dt*np.array(range(startMSD, endMSD)), MSD[startMSD:endMSD], 1)[3]
    # D = coeff[1]/(2*d*1000000)  # slope of MSD plot equals 2 D. Divide by 1 million to go from nm^2/ps to m^2/s
    times = dt*np.array(range(0,endMSD))
    avg_D, std_D = d_error(nMC, startMSD, endMSD, nT, limits, times, MSD, d)
    return MSD, endMSD, limits, avg_D, std_D


if __name__ == '__main__':
    args = initialize()
    # pos, _, trj_times, _ = gp('%s' % args.input, '%s' % args.comp, '%s' % args.LC_type, '%s' % args.solv)
    # f = open('d_pos', 'w')
    # np.save(f, pos)
    # f.close()
    # f = open('d_traj', 'w')
    # np.save(f, trj_times)
    # f.close()
    trj_times = np.load('d_traj')
    pos = np.load('d_pos')
    print 'got positions'
    nT = len(trj_times)
    Nbootstraps = args.nboot
    frontfrac = args.frontfrac
    fracshow = args.fracshow
    d = args.dim
    dt = (trj_times[len(trj_times) - 1] - trj_times[0])/(len(trj_times) - 1)
    print 'Getting the D ...'
    MSD, endMSD, limits, D_av, D_std = dconst(pos, nT, Nbootstraps, frontfrac, fracshow, d, dt, args.nMC)
    errorevery = int(np.ceil(fracshow*nT/100.0))  # plot only 100 bars total
    plt.errorbar(dt*np.array(range(0,endMSD)),MSD[:endMSD],yerr=[limits[0,:endMSD],limits[1,:endMSD]],errorevery=errorevery)
    plt.ylabel('MSD (nm^2)')
    plt.xlabel('time (ps)')
    plt.title('D = %s $\pm$ %s m$^2$/s' % (D_av, D_std))
    plt.show()