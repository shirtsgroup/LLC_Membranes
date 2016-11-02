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
    parser.add_argument('-b', '--nboot', default=2000, help='Number of bootstrap trials to be run')
    parser.add_argument('-f', '--frontfrac', default=0, help='Where to start fitting line for diffusivity calc')
    parser.add_argument('-F', '--fracshow', default=0.3, help='Percent of graph to show, also where to stop fitting line'
                                                             'during diffusivity calculation')
    parser.add_argument('-d', '--dim', default=1, help='Number of dimensions in which diffusion is being calculated')
    parser.add_argument('-l', '--LC_type', default='HII', help='Type of liquid crystal being studied')
    parser.add_argument('-s', '--solv', default='no', help='Is the system solvated or not?')
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


def dconst(x, nT, Nbootstraps, frontfrac, fracshow, d, dt):
    no_comp = len(x[0, :, 0])
    MSD, MSDs = msd(x, nT, no_comp)
    limits = bootstrap(nT, Nbootstraps, no_comp, MSDs, MSD)
    endMSD = int(np.floor(nT*fracshow))
    startMSD = int(np.floor(nT*frontfrac))
    coeff = Poly_fit.poly_fit(dt*np.array(range(startMSD, endMSD)), MSD[startMSD:endMSD], 1)[3]
    D = coeff[1]/(2*d*100)  # slope of MSD plot equals 2 D. Divide by 0.01 to go from nm^2/ps to cm^2/s
    return D, MSD, endMSD, limits


if __name__ == '__main__':
    args = initialize()
    pos, _, trj_times, _ = gp('%s' % args.input, '%s' % args.comp, '%s' % args.LC_type, '%s' % args.solv)
    nT = len(trj_times)
    Nbootstraps = args.nboot
    frontfrac = args.frontfrac
    fracshow = args.fracshow
    d = args.dim
    dt = (trj_times[len(trj_times) - 1] - trj_times[0])/(len(trj_times) - 1)
    D, MSD, endMSD, limits = dconst(pos, nT, Nbootstraps, frontfrac, fracshow, d, dt)
    errorevery = int(np.ceil(fracshow*nT/100.0))  # plot only 100 bars total
    plt.errorbar(dt*np.array(range(0,endMSD)),MSD[:endMSD],yerr=[limits[0,:endMSD],limits[1,:endMSD]],errorevery=errorevery)
    plt.ylabel('MSD')
    plt.xlabel('time (ps)')
    plt.title('Stackoverflow/Shirts algorithm MSD results')
    plt.show()