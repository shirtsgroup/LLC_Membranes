#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from past.utils import old_div
import argparse
import numpy as np
import matplotlib.pyplot as plt
import Poly_fit
import mdtraj as md
import top
import time

def initialize():

    parser = argparse.ArgumentParser(description='Calculate diffusion coefficient')
    parser.add_argument('-t', '--trajectory', default='wiggle.trr', help='Path to input file')
    parser.add_argument('-g', '--gro', default='wiggle.gro', help='Name of .gro coordinate file')
    parser.add_argument('-r', '--residue', default='SOL', help='Name of residue whose diffusivity we want')
    parser.add_argument('--itp', default='/usr/local/gromacs/share/gromacs/top/amber99.ff/tip3p.itp', help='Name of itp'
                        'describing topology of residue')
    parser.add_argument('-b', '--nboot', default=200, help='Number of bootstrap trials to be run')
    parser.add_argument('-f', '--frontfrac', default=0.05, type=float, help='Where to start fitting line on msd curve')
    parser.add_argument('-F', '--fracshow', default=.2, type=float, help='Percent of graph to show, also where to stop '
                        'fitting line during diffusivity calculation')
    parser.add_argument('-a', '--axis', default='xyz', type=str, help='Which axis to compute msd along')

    args = parser.parse_args()

    return args


def autocorrFFT(x):

    N = len(x)
    F = np.fft.fft(x, n=2*N)  # 2*N because of zero-padding
    PSD = F * F.conjugate()
    res = np.fft.ifft(PSD)
    res = (res[:N]).real   # now we have the autocorrelation in convention B
    n= N*np.ones(N)-np.arange(0, N)  # divide res(m) by (N-m)
    return old_div(res,n)  # this is the autocorrelation in convention A


def msd_fft(r):

    N = len(r)
    # N = r.shape[1]
    # print N
    D = np.square(r).sum(axis=1)
    D = np.append(D, 0)
    S2 = sum([autocorrFFT(r[:, i]) for i in range(r.shape[1])])
    Q = 2*D.sum()
    S1 = np.zeros(N)
    for m in range(N):
      Q = Q-D[m-1]-D[N-m]
      S1[m] = old_div(Q,(N-m))
    return S1-2*S2


def msd2(x):

    nT = x.shape[0]
    no_comp = x.shape[1]
    MSD = np.zeros([nT], dtype=float)
    MSDs = np.zeros([nT, no_comp], dtype=float)  # a set of MSDs per particle

    for n in range(no_comp):
        MSDs[:, n] = msd_fft(x[:, n, :])
        MSD += MSDs[:, n]
    MSD /= no_comp

    return MSD, MSDs


def msd(x, ndx):
    """
    Straightforward way to calculte msd
    :param x: positions of centers of mass of all particles for each frame, numpy array [nframes, natoms, dim]
    :param ndx: list of indices to include in msd calculation (x = 0, y = 1, z = 2)
    :return: Average MSD and individual particle MSDs
    """
    n = x.shape[1]  # number of atoms
    N = x.shape[0]  # number of frames

    MSD = np.zeros([N])
    MSDs = np.zeros([N, n])

    for m in range(N):  # there nT different length intervals we can look at
        for k in range(N - m - 1):  # there are N - m - 1 independent intervals of length m
            MSDs[m, :] += np.linalg.norm(x[k + m, :, ndx] - x[k, :, ndx], axis=0)**2  # mean square displacement of all particles summed over all intervals of length m
        MSDs[m, :] /= (N - m)  # divide by the number of intervals to get an average msd over length m
        MSD[m] = np.mean(MSDs[m, :])

    return MSD, MSDs


def bootstrap(nT, Nbootstraps, no_comp, MSDs, MSD, d, times, start, end):

    eMSDs = np.zeros([nT, Nbootstraps], dtype=float)  # a set of MSDs per particle

    # Now, bootstrap over number of particles, assuming that the particles are sufficiently independent
    for b in range(Nbootstraps):
        indices = np.random.randint(0, no_comp, no_comp)  # randomly select N of the particles with replacement
        # ^ makes a list of length Nparticles each with a random number from 0 to  Nparticles
        for n in range(no_comp):
            eMSDs[:, b] += MSDs[:, indices[n]]  # for this bootstrap trial b, add the MSDs of a randomly selected particle
            # to the particle at each time step
        eMSDs[:, b] /= no_comp  # Divide every timestep by Nparticles -- average the MSDs

    limits = np.zeros([2, nT], dtype=float)  # a set of MSDs per particle
    # now, let's determine a 95\% error bound for each tau (out of the
    # Nbootstrapped MSD's, use that for the error bars
    for t in range(0, nT):
        limits[0, t] = np.abs(np.percentile(eMSDs[t, :], 2.5) - MSD[t])
        limits[1, t] = np.abs(np.percentile(eMSDs[t, :], 97.5) - MSD[t])

    # generate error in slope
    tot_pts = int(end - start)
    # Weight least squares -- first generate a weight matrix
    W = np.zeros((tot_pts, tot_pts))
    for i in range(tot_pts):
        W[i, i] = old_div(1,((limits[0, i + start])**2))

    slopes = []
    for i in range(Nbootstraps):
        y_fit, _, slope_error, _, A = Poly_fit.poly_fit(times[start:end], eMSDs[start:end, i], 1, W)
        slopes.append(A[1])

    y_fit = Poly_fit.poly_fit(times[start:end], MSD[start:end], 1, W)[0]
    plt.plot(times[start:end], y_fit, '--', color='black', label='Linear Fit')

    return limits, np.mean(slopes)/(2*d*1000000), np.std(slopes)/(2*d*1000000)


def d_error(startfit, endMSD, nT, limits, times, MSD, d):

    tot_pts = (int(endMSD) - int(startfit))

    # Weight least squares -- first generate a weight matrix
    W = np.zeros((tot_pts, tot_pts))
    for i in range(tot_pts):
        W[i, i] = old_div(1,((limits[0, i + startfit])**2))

    y_fit, _, slope_error, _, A = Poly_fit.poly_fit(times[startfit:endMSD], MSD[startfit:endMSD], 1, W)
    plt.plot(times[startfit:endMSD], y_fit, '--', color='black', label='Linear Fit')

    return A[1]/(2*d*1000000), slope_error/(2*d*1000000)


def dconst(x, nT, Nbootstraps, frontfrac, fracshow, d, dt, ndx):

    start = time.time()
    MSD, MSDs = msd(x, ndx)
    end = time.time()
    print('MSDs done in %1.2f seconds' % (end - start))

    endMSD = int(nT*fracshow)
    startMSD = int(nT*frontfrac)

    times = dt*np.array(list(range(0, endMSD)))

    limits, avg_D, std_D = bootstrap(nT, Nbootstraps, x.shape[1], MSDs, MSD, d, times, startMSD, endMSD)
    print('Bootstrapping Done')

    # avg_D, std_D = d_error(startMSD, endMSD, nT, limits, times, MSD, d)

    return MSD, endMSD, limits, avg_D, std_D


if __name__ == '__main__':

    args = initialize()  # initialize variables passed to the script

    t = md.load(args.trajectory, top=args.gro)  # load trajectory
    nT = t.n_frames  # number of frames

    res = args.residue
    if res == 'SOL':  # mdtraj changes the name from SOL to HOH
        res = 'HOH'

    selection = [a.index for a in t.topology.atoms if a.residue.name == res]

    pos = t.xyz[:, selection, :]  # positions of all atoms of interest

    topology = top.Top(args.itp)  # read topology
    atoms_per_residue = topology.natoms  # number atoms in a single residue
    matoms = topology.atom_masses  # mass of the atoms in residue
    mres = np.sum(matoms)  # total mass of residue

    dimension = len(args.axis)  # number of dimensions in which msd is being computed
    ndx = []
    if 'x' in args.axis:
        ndx.append(0)
    if 'y' in args.axis:
        ndx.append(1)
    if 'z' in args.axis:
        ndx.append(2)

    com = np.zeros([nT, pos.shape[1]//atoms_per_residue, 3])  # track the center of mass of each residue

    for f in range(nT):
        for i in range(com.shape[1]):
            w = (pos[f, i*atoms_per_residue:(i+1)*atoms_per_residue, :].T * matoms).T  # weight each atom in the residue by its mass
            com[f, i, :] = np.sum(w, axis=0) / mres  # sum the coordinates and divide by the mass of the residue

    dt = t.time[-1] - t.time[-2]  # time step (assuming equispaced time points)

    MSD, MSDs = msd(pos, ndx)
    startMSD = int(nT*args.frontfrac)
    endMSD = int(nT*args.fracshow)
    limits, D_avg, D_std = bootstrap(nT, args.nboot, pos.shape[1], MSDs, MSD, len(ndx), t.time, startMSD, endMSD)
    # MSD, endMSD, limits, D_av, D_std = dconst(com, nT, args.nboot, args.frontfrac, args.fracshow, dimension, dt, ndx)
    # errorevery = int(np.ceil(args.fracshow*nT/100.0))  # plot only 100 bars total
    errorevery = int(np.ceil(nT/100.0))  # plot only 100 bars total
    # plt.errorbar(dt*np.array(list(range(0, endMSD))), MSD[:endMSD], yerr=[limits[0, :endMSD], limits[1, :endMSD]],
    #              errorevery=errorevery, label='Calculated MSD')
    plt.errorbar(dt*np.array(list(range(0, nT-1))), MSD[:-1], yerr=[limits[0, :-1], limits[1, :-1]],
             errorevery=errorevery, label='Calculated MSD')
    plt.ylabel('MSD ($nm^2$)', fontsize=14)
    plt.xlabel('time (ps)', fontsize=14)
    plt.gcf().get_axes()[0].tick_params(labelsize=14)
    plt.title('D = %1.2e $\pm$ %1.2e $m^{2}/s$' % (D_avg, D_std))
    plt.legend(loc=2)
    plt.tight_layout()
    plt.savefig('Diffusivity_%s.png' % args.axis)
    plt.show()
