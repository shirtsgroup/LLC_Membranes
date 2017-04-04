#! /usr/bin/env python

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

import argparse
from Get_Positions import get_positions
import Atom_props
import matplotlib.pyplot as plt
import time
import Poly_fit
import numpy as np
import Diffusivity
import Concentration
import subprocess
import mdtraj as md
from llclib import physical


def initialize():

    parser = argparse.ArgumentParser(description = 'Caclulate Ionic Conductivity using the Nernst Einstein relation or'
                                                   'the Collective Diffusion Model')

    parser.add_argument('-t', '--traj', default='wiggle.trr', type=str, help='Trajectory file (.xtc or .trr)')
    parser.add_argument('-g', '--gro', default='wiggle.gro', type=str, help='Coordinate file')
    parser.add_argument('-d', '--build_mon', default='NAcarb11V', type=str, help='Monomer with which the system was built')
    parser.add_argument('-i', '--ion', default='NA', help='Name of ion(s) being used to calculate ionic conductivity')
    parser.add_argument('-b', '--buffer', default='.2', help='Distance into membrane from min and max where current '
                                                            'measurements will be made')
    parser.add_argument('-T', '--temp', default=300, help='System Temperature, Kelvin')
    parser.add_argument('-B', '--nboot', default=200, help='Number of bootstrap trials to be run')
    parser.add_argument('-m', '--nMC', default=1000, help='Number of Monte Carlo trials to estimate error in D and Dq')
    parser.add_argument('-a', '--arrays', default='off', help='If positions, id array are already saved')
    parser.add_argument('-S', '--suffix', default='saved', help='Suffix to append to position and id arrays when saving')
    parser.add_argument('-r', '--frontfrac', default=0.05, help='Where to start fitting line for diffusivity calc')
    parser.add_argument('-F', '--fracshow', default=.8, help='Percent of graph to show, also where to stop fitting line'
                                                             'during diffusivity calculation')
    #parser.add_argument('-t', '--nTsub', default=20, help='Number of subtrajectories to break into for generating stats')
    parser.add_argument('-M', '--method', default='B', help='Which method to use to calculate Ionic Conductivity: CD ='
                                                             'Collective Diffusion, NE = Nernst Einstein, B = both')
    parser.add_argument('-l', '--load', help='Load arrays if they were saved from a previous calculation',
                        action="store_true")
    parser.add_argument('-s', '--save', help='Save arrays generated for future use', action="store_true")

    args = parser.parse_args()

    return args


def nernst_einstein(D, D_std, C, C_std, T):
    """

    :param D:
    :param D_std:
    :param C:
    :param C_std:
    :param T:
    :return: Ionic Conductivity as calculated by the nernst einsten relation
    """
    NE_av = q**2*C*D/(kb*float(T))
    NE_error = NE_av*np.sqrt((D_std/D)**2 + (C_std/C)**2)

    return NE_av, NE_error

if __name__ == '__main__':

    start = time.time()
    args = initialize()

    if args.load:

        # If we are just looking at different trajectories or fitting different parts of msd curves, we can just load what
        # has already been calculated

        pos = np.load('pos_%s' % args.suffix)

        id = np.load('identity_%s' % args.suffix)

        dt = np.load('dt_%s' % args.suffix)

        Conc_props = np.load('C_props_%s' % args.suffix)
        C = Conc_props[0]
        C_std = Conc_props[1]
        Lc = Conc_props[2]
        conv = Conc_props[3]
        z_2 = Conc_props[4]
        z_1 = Conc_props[5]

        D_tails = np.load('D_tails_%s' % args.suffix)
        D_av = D_tails[0]
        D_std = D_tails[1]

        no_comp = np.shape(pos)[1]
        nT = np.shape(pos)[2]

        print 'All arrays/properties loaded'

    else:

        t = md.load(args.traj, top=args.gro)
        keep = [a.index for a in t.topology.atoms if a.name == args.ion]
        t_ion = t.atom_slice(keep)
        pos = t.xyz
        print pos.shape
        pos_ion = t_ion.xyz
        time = t.time

        C, C_std, cross, Lc, z_2, z_1 = physical.conc(t, t_ion, args.buffer)

        factor = 1*10**9  # number of nm in a m
        conv = (Lc/factor)/(cross/(factor**2))

        Conc_props = np.array([C, C_std, Lc, conv, z_1, z_2])

        no_comp = np.shape(pos)[1]
        nT = np.shape(pos)[2]

        dt = time[1] - time[0]

        id = np.array([a.name for a in t.topology.atoms])

        if args.save:

            with open('pos_%s' % args.suffix, 'w') as f: # Save the positions for running the script again
                np.save(f, pos)

            with open('%spos_%s' % (args.ion, args.suffix), 'w') as f:  # Save the positions for running the script again
                np.save(f, pos_ion)

            with open('identity_%s' % args.suffix, 'w') as f:
                np.save(f, id)

            with open('C_props_%s' % args.suffix, 'w') as f:
                np.save(f, Conc_props)

            with open('dt_%s' % args.suffix, 'w') as f:
                np.save(f, dt)

            print 'All arrays/properties saved'

# Constants

kb = 1.38*10**-23  # Boltzmann constant [=] J/K
e = 1.602 * 10 ** -19  # elementary charge (C)
q = Atom_props.charge[args.ion] * e  # elementary charge (1.602e-19 C) * charge of ion
frontfrac = args.frontfrac
fracshow = args.fracshow
d = 3  # dimensionality


if args.method == 'NE' or args.method == 'B':

    if args.arrays == 'off':

        print 'Calculating Diffusivity'
        # Calculate diffusivity (m^2/s)
        _, _, _, D_av, D_std = Diffusivity.dconst(pos_ion, nT, args.nboot, frontfrac, fracshow, d, dt, args.nMC)
        D_tails = [D_av, D_std]  # Get it? D_tails, i.e. details about D

        if args.save:
            with open('D_tails_%s' % args.suffix, 'w') as f:
                np.save(f, D_tails)

    IC_NE, IC_NE_std = nernst_einstein(D_av, D_std, C, C_std, '%s' % args.temp)
    print 'Nernst Einstein Ionic Conductivity: %s +/- %s S/m' % (IC_NE, IC_NE_std)

    if args.method == 'NE':
        exit()


def dQ(frame, positions, channel_length, zmax, zmin):
    q = 0
    for i in range(np.shape(positions)[1]):  # find out how far each ion has moved in the z direction this frame
        curr = positions[2, i, frame]  # current z position
        prev = positions[2, i, frame - 1]  # previous z position
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
        if prev_frame == zmin and current_frame == zmin:
            pass
        if prev_frame == zmax and current_frame == zmax:
            pass
        displacement = current_frame - prev_frame
        # q += e*(Atom_props.charge[id[0, i, frame]])*displacement/channel_length
        q += e*(Atom_props.charge[id[0, i]])*displacement/channel_length
    return q

nT -= 1
n_sub = args.nTsub  # number of sub-trajectories to used for error analysis
nT_sub = nT/n_sub  # number of points in that subinterval

if args.arrays == 'off':

    print 'Calculating q'
    dq_all = np.zeros([nT])
    for i in range(nT):
        dq_all[i] = dQ(i, pos, Lc, z_2, z_1)

    f = open('dq_%s' % args.suffix, 'w')
    np.save(f, dq_all)
    f.close()
    print 'q saved'

else:
    dq_all = np.load('dq_%s' % args.suffix)
    print 'q loaded'

# Break dq_all into sub-trajectories
dq = np.zeros([n_sub, nT_sub])
for i in range(n_sub):
    dq[i, :] = dq_all[i * nT_sub: (i + 1) * nT_sub]

dq_cum = np.zeros([n_sub, nT_sub])
dq_cum[:, 0] = dq[:, 0]
for i in range(nT_sub):
    dq_cum[:, i] = dq[:, i] + dq_cum[:, i - 1]

msd = np.zeros([n_sub, nT_sub])
itau = 0  # counts up to the number of trajectory frames
while itau < nT_sub:
    ncount = (nT_sub-itau)
    for t in range(nT_sub-itau):  # run for total number of time points from itau to nT
        for i in range(n_sub):
            xo = dq_cum[i, t + itau] - dq_cum[i, t]
            msd[i, itau] += np.dot(xo, xo)  # Should there be a division by 't' in here?
    msd[:, itau] /= ncount
    itau += 1

avg_std = np.zeros([2, nT_sub])
for i in range(nT_sub):
    avg_std[0, i] = np.mean(msd[:, i])
    avg_std[1, i] = np.std(msd[:, i])

times = np.linspace(dt, nT_sub*dt, nT_sub)
startfit = int(np.floor(nT_sub*frontfrac))
endMSD = int(np.floor(nT_sub*fracshow))

tot_pts = endMSD - startfit
slopes = np.zeros([n_sub])
for i in range(n_sub):
    plt.plot(times[:endMSD], msd[i, :endMSD])
    y_fit, _, slope_error, _, A = Poly_fit.poly_fit(times[startfit:endMSD], msd[i, startfit:endMSD], 1, 'none')
    # plt.plot(times[startfit:endMSD], y_fit)
    slopes[i] = A[1]

plt.figure()
plt.title('%s MSD trajectories' % n_sub)
plt.xlabel('Time (ps)')
plt.ylabel('MSD (e$^2$)')

ic_cd_avg = (np.mean(slopes)*(1*10**12)*conv)/(2*kb*float(args.temp))
ic_cd_std = (np.std(slopes)*(1*10**12)*conv)/(2*kb*float(args.temp))
ic_cds = slopes*((1*10**12)*conv/(2*kb*float(args.temp)))

# Bootstrap the slopes
from random import randint

nboot = 200000
means = np.zeros([nboot])
for i in range(nboot):
    sum = 0
    for j in range(n_sub):
        index = randint(0, n_sub - 1)
        sum += slopes[index]
    means[i] = sum / n_sub
means *= ((1*10**12)*conv)/(2*kb*float(args.temp))
bootstrap_mean = np.mean(means)
bootstrap_std = np.std(means)

plt.figure()
plt.ylabel('Count')
plt.xlabel('')
plt.hist(means, bins=100)
print "Ionic Conductivity: %s +/- %s S/m" % (ic_cd_avg, ic_cd_std)
print "Bootstrapped Ionic Conductivity: %s +/- %s S/m" % (bootstrap_mean, bootstrap_std)


plt.show()
# plt.errorbar(times[:endMSD], avg_std[0, :endMSD], yerr=avg_std[1, :endMSD])
# plt.show()
exit()
# Monte Carlo Method:
# Fit a line through the points
# y_fit, r_squared, std2, coeff = Poly_fit.poly_fit(times[startfit:endMSD], avg_std[0, startfit:endMSD], 1)
# Generate uncertainty in the slope of the plot by assuming a gaussian distribution in the error bars. Randomly select
# a number centered at the plot point from avg_std[0, :] with standard deviation defined entries in avg_std[1, :], then
# take the avg/std of this for the final error calculation
# nMC = 1000  # number of Monte Carlo trajectories
# pts = np.zeros([nMC, (endMSD - startfit)])
# for i in range(nMC):
#     for j in range(0, (endMSD - startfit)):
#         pts[i, j] = avg_std[1, j + startfit]*np.random.randn() + avg_std[0, j + startfit]
#
# slopes = np.zeros([nMC])
# for i in range(nMC):
#     slopes[i] = Poly_fit.poly_fit(times[startfit:endMSD], pts[i, :], 1)[3][1]
#
# avg_slope = np.mean(slopes)
# std_slope = np.std(slopes)


# Weight least squares -- first generate a weight matrix
tot_pts = endMSD - startfit
W = np.zeros((tot_pts, tot_pts))
for i in range(tot_pts):
    W[i, i] = 1/((avg_std[1, i + startfit])**2)
print np.shape(W)
y_fit, _, slope_error, _, A = Poly_fit.poly_fit(times[startfit:endMSD], avg_std[0, startfit:endMSD], 1, W)

avg_slope = A[1]

ic_cd_avg = (avg_slope*(1*10**12)*conv)/(2*kb*float(args.temp))
ic_cd_std = (slope_error*(1*10**12)*conv)/(2*kb*float(args.temp))
print "Ionic Conductivity: %s +/- %s S/m" % (ic_cd_avg, ic_cd_std)
plt.suptitle('MSD of Q vs time')
plt.title('Ionic Conductivity: %s +/- %s S/m' % (ic_cd_avg, ic_cd_std))
plt.xlabel('Time (ps)')
plt.ylabel('MSD (e$^2$)')
plt.plot(times[:endMSD], avg_std[0, :endMSD])
plt.plot(times[startfit:endMSD], y_fit)
end = time.time()
print 'total execution time: %s seconds' % (end - start)
plt.show()