#! /usr/bin/env python

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

from builtins import range
from past.utils import old_div
import argparse
import Atom_props
import matplotlib.pyplot as plt
import time
import Poly_fit
import numpy as np
import Diffusivity
import mdtraj as md
from llclib import physical
import tqdm

"""
The purpose of this script is to calculate the ionic conductivity of a given LLC membrane using the Nernst Einstein
relation (Method 1) and the Collective Diffusion model (Method 2), the latter being described in the following paper:

    Y. Liu and F. Zhu "Collective diffusion model for ion conduction through microscopic channels," Biophysical Journal,
    vol 104, no. 2, pp. 368-376, Jan. 2013.

Ionic conductivity quantifies the movement of ions from site to site through defects in the crystal lattice of
ionic solids or aqueous solutions.

METHOD 1: Nernst Einstein

METHOD 2: Collective Diffusion

To calculate the ionic conductivity, we need to know the current of ions. For current through a channel we define the
current through the channel as:

    I = sum ( e_i * v_i / Lc )                                                                                (1)

where e and v are the charge and velocity of components i located in the channel. Lc is the channel length defined by
two planes at z1 and z2. Thus Lc = z2 - z1. Equation (1), on average, is the same as the current through the entire
system.

The charge transfer through the channel is defined as:

    Q(t) = integral[0, t](I(t)dt)                                                                             (2)

One can use the equilibrium relation between the mean square displacement of Q and its velocity autocorrelation function

    <Q^2(t)> = 2*integral[0, t]dt'(t - t')<I(0)I(t')>                                                         (3)

The equilibrium diffusion coefficient of Q, Dq, is given by:

    Dq = integral[0, infinity](dt*expected_value(I(0)*I(t)))                                                  (4)

Alternatively, we can plot Q(t), find its slope in the linear region and get Dq (much like calculating diffusion
coefficients). Define the change in Q between each time step as:

    delta_Q = sum(e_i*delta_z/Lc)                                                                             (5)

Where we sum over atoms located within the channel region bounded by z_max and z_min. e_i is the charge of atom i.
delta_z is the net displacement of an atom. If an atom exits or enters the boundaries between time steps, delta_z is
taken only as the displacement of the atom within the channel region. For LLC membranes, the channel membrane is taken
as the entire membrane cross section bounded by z_max and z_min since we are looking at bulk properties as opposed to
pure channel conduction. We then cumulate the delta_Qs to construct Q(t). Finding the slope of the linear region will
give the value for Dq. Simulations need to be run for a long time to get accurate values. (LLC systems should be run for
a microsecond or longer to give accurate values)

The steady state current is given by:

    I_steady = (Dq * V) / (k_b * T)                                                                           (6)

Solving for I/V gives the system conductance:

    c = Dq / (k_b * T)                                                                                        (7)

This can be used to calculate the ionic conductance with an equilibrium simulation

ALTERNATIVE METHODS:

We can calculate conductance directly from current and voltage if an electric field is applied across the
membrane. To do so, calculate the average current using a time average of equation (1) and the potential difference
across the membrane (which can be found using gmx potential with a careful look at the output .xvg file).

One can also employ the Computational Electrophysiology Module included with GROMACS to create a double layer membrane
system which maintains a concentration gradient across the membrane, inducing a current.
"""


def initialize():

    parser = argparse.ArgumentParser(description='Calculate Ionic Conductivity using the Nernst Einstein relation or'
                                                 'the Collective Diffusion Model')

    parser.add_argument('-t', '--traj', default='wiggle.trr', type=str, help='Trajectory file (.xtc or .trr)'
                                                                             'IMPORTANT: Pre-process the trajectory with gmx trjconv -pbc nojump')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate file')
    parser.add_argument('-d', '--build_mon', default='NAcarb11V', type=str, help='Monomer with which the system was '
                                                                                 'built')
    parser.add_argument('-i', '--ion', default='NA', help='Name of ion(s) being used to calculate ionic conductivity')
    parser.add_argument('-b', '--buffer', default=0, type=float, help='Distance into membrane from min and max where '
                                                                      'current measurements will be made')
    parser.add_argument('-T', '--temp', default=300, help='System Temperature, Kelvin')
    parser.add_argument('-B', '--nboot', default=200, type=int, help='Number of bootstrap trials to be run')
    parser.add_argument('-m', '--nMC', default=1000, help='Number of Monte Carlo trials to estimate error in D and Dq')
    parser.add_argument('-S', '--suffix', default='saved',
                        help='Suffix to append to position and id arrays when saving')
    parser.add_argument('-r', '--frontfrac', default=0.16, type=float, help='Where to start fitting line for diffusivity'
                                                                           'calc')
    parser.add_argument('-F', '--fracshow', default=0.5, type=float, help='Percent of graph to show, also where to stop'
                                                                          'fitting line during diffusivity calculation')
    parser.add_argument('--nTsub', default=5, type=int, help='Number of subtrajectories to break into for generating '
                                                              'stats')
    parser.add_argument('-M', '--method', default='NE', help='Which method to use to calculate Ionic Conductivity: CD ='
                                                             'Collective Diffusion, NE = Nernst Einstein, B = both')
    parser.add_argument('-l', '--load', help='Load arrays if they were saved from a previous calculation',
                        action="store_true")
    parser.add_argument('--discard', type=int, help='Specify the number of nanoseconds to discard starting'
                                                    'from the beginning of the simulation')
    parser.add_argument('--noshow', action="store_true", help='Specify this to not show any plots')
    parser.add_argument('-a', '--axis', default='xyz', type=str, help='Which axis to compute msd along')

    args = parser.parse_args()

    return args


def nernst_einstein(D, D_std, C, C_std, T):
    """
    :param D: calculated diffusivity
    :param D_std: error in calculated diffusivity
    :param C: concentration of ions in membrane
    :param C_std: error in concentration
    :param T: temperature simulation is run at
    :return: Ionic Conductivity as calculated by the nernst einsten relation
    """

    NE_av = q ** 2 * C * D / (kb * float(T))
    NE_error = NE_av * np.sqrt((old_div(D_std, D)) ** 2 + (old_div(C_std, C)) ** 2)

    return NE_av, NE_error


def dQ(frame, positions, channel_length, zmax, zmin, charges, id):

    q = 0

    for i in range(positions.shape[1]):  # find out how far each ion has moved in the z direction this frame
        curr = positions[frame, i, 2]  # current z position
        prev = positions[frame - 1, i, 2]  # previous z position
        if zmin <= curr <= zmax:  # check that the ion is in the channel over which we are measuring
            current_frame = curr  # If it, that is its position at the current frame
        elif curr > zmax:  # If the particle is above zmax, we are still potentially interested in it
            current_frame = zmax
        elif curr < zmin:  # If the particle is below zmin, we are still potentially interested in it
            current_frame = zmin
        if zmin <= prev <= zmax:  # Do the same checks on the particle from the previous frame
            prev_frame = prev
        elif prev > zmax:
            prev_frame = zmax
        elif prev < zmin:
            prev_frame = zmin
        displacement = current_frame - prev_frame
        q += e * (charges[id[i]]) * displacement / channel_length

    return q


def dQ2(z, zmax, zmin, charges, id):
    """
    Calculate the displacement of all positions from frame to frame
    :param z: trajectory of positions [nframes, npoints, 1]  --> z-direction only
    :return: a trajectory n(t) describing the delta n at each time step
    """
    q = np.zeros(z.shape[0])

    for t in tqdm.tqdm(range(1, q.shape[0])):
        for i in range(z.shape[1]):
            current = z[t, i]
            previous = z[t - 1, i]
            if zmax >= current >= zmin:
                if zmax >= previous >= zmin:
                    q[t] += (current - previous)*charges[id[i]]
                elif previous > zmax:
                    q[t] += (current - zmax)*charges[id[i]]
                elif previous < zmin:
                    q[t] += (current - zmin)*charges[id[i]]
            elif current > zmax:
                if zmax >= previous >= zmin:
                    q[t] += (zmax - previous)*charges[id[i]]
            elif current < zmin:
                if zmax >= previous >= zmin:
                    q[t] += (zmin - previous)*charges[id[i]]

    q *= (e / (zmax - zmin))

    return q


if __name__ == '__main__':

    start = time.time()
    args = initialize()

    # Constants

    kb = 1.38 * 10 ** -23  # Boltzmann constant [=] J/K
    e = 1.602 * 10 ** -19  # elementary charge (C)
    q = Atom_props.charge[args.ion] * e  # elementary charge (1.602e-19 C) * charge of ion
    frontfrac = args.frontfrac
    fracshow = args.fracshow

    d = len(args.axis)  # number of dimensions in which msd is being computed
    ndx = []
    if 'x' in args.axis:
        ndx.append(0)
    if 'y' in args.axis:
        ndx.append(1)
    if 'z' in args.axis:
        ndx.append(2)

    t = md.load(args.traj, top=args.gro)

    keep = [a.index for a in t.topology.atoms if a.name == args.ion]

    pos = t.xyz
    pos_ion = t.xyz[:, keep, :]
    time = t.time

    C, C_std, cross, Lc, z_2, z_1 = physical.conc(t, args.ion, args.buffer)

    factor = 1 * 10 ** 9  # number of nm in a m
    conv = old_div((old_div(Lc, factor)), (old_div(cross, (factor ** 2))))

    no_comp = pos.shape[1]  # np.shape(pos)[1]
    nT = t.n_frames  # np.shape(pos)[0]

    dt = t.timestep  # time[1] - time[0]

    id = np.array([a.name for a in t.topology.atoms])

    if args.method == 'NE' or args.method == 'B':

        # if not args.load:

        print('Calculating Diffusivity')

        MSD, endMSD, limits, D_avg, D_std = Diffusivity.dconst(pos_ion, nT, args.nboot, frontfrac, fracshow, d, dt, ndx)

        fit = 0
        while fit == 0:

            endMSD = int(nT*fracshow)
            startMSD = int(nT*frontfrac)
            D_avg, D_std = Diffusivity.d_error(startMSD, endMSD, nT, limits, t.time, MSD, d)
            errorevery = int(np.ceil(nT / 100.0))  # plot only 100 bars total
            plt.errorbar(dt * np.array(list(range(0, nT - 1))), MSD[:-1], yerr=[limits[0, :-1], limits[1, :-1]],
                         errorevery=errorevery, label='Calculated MSD')
            plt.ylabel('MSD ($nm^2$)', fontsize=14)
            plt.xlabel('time (ps)', fontsize=14)
            plt.gcf().get_axes()[0].tick_params(labelsize=14)
            plt.title('D = %1.2e $\pm$ %1.2e $m^{2}/s$' % (D_avg, D_std))
            plt.legend(loc=2)
            plt.tight_layout()
            plt.savefig('Diffusivity_%s.png' % args.axis)
            plt.ion()
            plt.show()
            fit = int(input("Type '1' if the fit looks good: "))
            if fit != 1:
                print('Press enter to following prompts to leave as is')
                frontfrac = float(input("Fraction into simulation to start fit: ") or frontfrac)
                fracshow = float(input("Fraction into simulation to stop fit: ") or fracshow)
            plt.clf()

        IC_NE, IC_NE_std = nernst_einstein(D_avg, D_std, C, C_std, '%s' % args.temp)
        print('Nernst Einstein Ionic Conductivity: {:1.2e} +/- {:1.2e} S/m'.format(IC_NE, IC_NE_std))

        if args.method == 'NE':
            exit()

    if not args.load:

        charge = Atom_props.charges(args.build_mon)
        charge['NA'] = 1.00

        print('Calculating q')
        # dq_all = np.zeros([nT - 1])
        # for i in tqdm.tqdm(range(nT - 1)):
        #     dq_all[i] = dQ(i + 1, pos, Lc, z_2, z_1, charge, id)
        dq_all = dQ2(pos[:, :, 2], z_2, z_1, charge, id)

        np.save('dq_%s' % args.suffix, dq_all)
        print('q saved')

    else:

        dq_all = np.load('dq_%s.npy' % args.suffix)
        print('q loaded')

    if args.discard:
        discarded_frames = int(
            old_div(args.discard, (old_div(dt, 1000))))  # convert nanoseconds to frames of discarded simulation
        print('First %s frames discarded (%s ns)' % (discarded_frames, args.discard))
        nT -= discarded_frames
    else:
        discarded_frames = 0

    n_sub = args.nTsub  # number of sub-trajectories to used for error analysis
    nT_sub = old_div((nT - 1), n_sub)  # number of points in that subinterval
    l_sub = dt * nT_sub / 1000  # ns in the subinterval

    # Break dq_all into sub-trajectories
    dq = np.zeros([n_sub, nT_sub])
    for i in range(n_sub):
        dq[i, :] = dq_all[i * nT_sub + discarded_frames: (i + 1) * nT_sub + discarded_frames]

    # cumulate q in each time interval
    dq_cum = np.zeros([n_sub, nT_sub])
    dq_cum[:, 0] = dq[:, 0]
    for i in range(nT_sub):
        dq_cum[:, i] = dq[:, i] + dq_cum[:, i - 1]

    msd = np.zeros([n_sub, nT_sub])
    itau = 0  # counts up to the number of trajectory frames
    while itau < nT_sub:
        ncount = (nT_sub - itau)
        for t in range(nT_sub - itau):  # run for total number of time points from itau to nT
            for i in range(n_sub):
                xo = dq_cum[i, t + itau] - dq_cum[i, t]
                msd[i, itau] += np.dot(xo, xo)  # Should there be a division by 't' in here?
        msd[:, itau] /= ncount
        itau += 1

    avg_std = np.zeros([2, nT_sub])
    for i in range(nT_sub):
        avg_std[0, i] = np.mean(msd[:, i])
        avg_std[1, i] = np.std(msd[:, i])

    times = np.linspace(0, (nT_sub - 1) * dt, nT_sub)

    plt.figure(1)

    avg = np.mean(msd, axis=0)
    std = np.std(msd, axis=0)

    startfit = int(np.floor(nT_sub * frontfrac))
    endMSD = int(np.floor(nT_sub * fracshow))
    start = startfit*dt/1000
    end = endMSD*dt/1000

    fit = 0
    while fit == 0:

        y_fit = Poly_fit.poly_fit(times[startfit:endMSD], avg[startfit:endMSD], 1, 'none')[0]
        plt.plot(times[startfit:endMSD], y_fit, '--', color='black')
        plt.errorbar(times[:nT_sub], avg, yerr=std, label='Calculated MSD')
        plt.xlabel('Time (ps)', fontsize=14)
        plt.ylabel('MSD (e$^2$)', fontsize=14)
        plt.gcf().get_axes()[0].tick_params(labelsize=14)
        plt.title('%s MSD trajectories' % n_sub)
        plt.legend(loc=2)
        plt.tight_layout()
        plt.savefig('MSD_CD.png')
        plt.ion()
        plt.show()
        fit = int(input("Type '1' if the fit looks good: "))
        if fit != 1:
            print('Press enter to following prompts to leave as is')
            start = float(input("Time (ns) into simulation to start fit: ") or start)
            startfit = int(1000*start/dt)
            end = float(input("Time (ns) into simulation to stop fit: ") or end)
            endMSD = int(1000*end/dt)
            plt.cla()

    # fit slopes to all sub-trajectories for error estimation
    tot_pts = endMSD - startfit
    slopes = np.zeros([n_sub])
    for i in range(n_sub):
        y_fit, _, slope_error, _, A = Poly_fit.poly_fit(times[startfit:endMSD], msd[i, startfit:endMSD], 1, 'none')
        slopes[i] = A[1]

    ic_cd_avg = old_div((np.mean(slopes) * (1 * 10 ** 12) * conv), (2 * kb * float(args.temp)))
    ic_cd_std = old_div((np.std(slopes) * (1 * 10 ** 12) * conv), (2 * kb * float(args.temp)))
    ic_cds = slopes * ((1 * 10 ** 12) * conv / (2 * kb * float(args.temp)))

    # Bootstrap the slopes
    from random import randint

    nboot = 200000
    means = np.zeros([nboot])
    for i in range(nboot):
        sum = 0
        for j in range(n_sub):
            index = randint(0, n_sub - 1)
            sum += slopes[index]
        means[i] = old_div(sum, n_sub)
    means *= old_div(((1 * 10 ** 12) * conv), (2 * kb * float(args.temp)))
    bootstrap_mean = np.mean(means)
    bootstrap_std = np.std(means)

    plt.figure(2)
    plt.ylabel('Count')
    plt.xlabel('')
    plt.hist(means, bins=100)

    print("Slope fit from %1.0f ns to %1.0f ns of each subtrajectory" % (start, end))
    print('Bootstrapped Ionic Conductivity: {:1.2e} +/- {:1.2e} S/m'.format(bootstrap_mean, bootstrap_std))

    if not args.noshow:
        plt.show()