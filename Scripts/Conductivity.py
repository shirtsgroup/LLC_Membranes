#!/usr/bin/python

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

parser = argparse.ArgumentParser(description = 'Run Cylindricity script')

parser.add_argument('-f', '--trajfile', default='wiggle_traj_754.gro', help = 'Path to input file')
parser.add_argument('-i', '--ion', default='NA', help='Name of ion(s) being used to calculate ionic conductivity')
parser.add_argument('-l', '--LC_type', default='HII', help='Type of liquid crystal. Should match that defined in '
                                                           'LC_class.py')
parser.add_argument('-s', '--solv', default='no', help='Is the system solvate?')
parser.add_argument('-b', '--buffer', default='.2', help='Distance into membrane from min and max where current '
                                                        'measurements will be made')
parser.add_argument('-T', '--temp', default=300, help='System Temperature, Kelvin')
parser.add_argument('-B', '--nboot', default=200, help='Number of bootstrap trials to be run')
parser.add_argument('-m', '--nMC', default=1000, help='Number of Monte Carlo trials to estimate error in D and Dq')
parser.add_argument('-a', '--arrays', default='on', help='If positions, id array are already saved')
parser.add_argument('-S', '--suffix', default='array612ns', help='Suffix to append to position and id arrays when saving')
parser.add_argument('-r', '--frontfrac', default=0.05, help='Where to start fitting line for diffusivity calc')
parser.add_argument('-F', '--fracshow', default=.4, help='Percent of graph to show, also where to stop fitting line'
                                                         'during diffusivity calculation')
parser.add_argument('-t', '--nTsub', default=10, help='Number of subtrajectories to break into for generating stats')
parser.add_argument('-M', '--method', default='B', help='Which method to use to calculate Ionic Conductivity: CD ='
                                                         'Collective Diffusion, NE = Nernst Einstein, B = both')

args = parser.parse_args()
start = time.time()
if args.arrays == 'on':

    # If we are just looking at different trajectories or fitting different parts of msd curves, we can just load what
    # has already been calculated

    pos = np.load('pos_%s' % args.suffix)
    print 'Positions loaded'

    id = np.load('identity_%s' % args.suffix)
    print 'IDs loaded'

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

    print 'All other constants loaded'

else:

    print 'Calculating Concentration first for memory reasons'
    C_list = subprocess.check_output(["python", "Concentration.py", "-i", "%s" % args.trajfile, "-c", "%s" % args.ion,
                                      "-b", "%s" % args.buffer, "-l", "%s" % args.LC_type, "-s" "%s" % args.solv]).splitlines()
    C_list = np.array([float(i) for i in C_list])

    C = C_list[0]
    C_std = C_list[1]
    cross = C_list[2]
    z_2 = C_list[3]
    z_1 = C_list[4]
    Lc = z_2 - z_1
    factor = 1*10**9  # number of nm in a m
    conv = (Lc/factor)/(cross/(factor**2))

    print 'Concentration Calculated'

    Conc_props = np.array([C_list[0], C_list[1], Lc, conv, C_list[4], C_list[3]])

    f = open('C_props_%s' % args.suffix, 'w')
    np.save(f, Conc_props)
    f.close()

    print 'Concentration properties saved (average, std, cross section, thickness, z boundaries)'

    pos, vel, trj_times, box = get_positions('%s' % args.trajfile, 'sys', '%s' % args.LC_type, '%s' % args.solv)

    f = open('pos_%s' % args.suffix, 'w')  # Save the positions for running the script again
    np.save(f, pos)
    f.close()
    print 'Atomic postions got and saved'

    pos_ion, vel, trj_times, box = get_positions('%s' % args.trajfile, '%s' % args.ion, '%s' % args.LC_type, '%s' % args.solv)

    f = open('%spos_%s' % (args.ion, args.suffix), 'w')  # Save the positions for running the script again
    np.save(f, pos_ion)
    f.close()
    print 'Ion positions got and saved'

    # Temporary Section for debugging
    # pos = np.load('pos_%s' % args.suffix)
    # print 'Positions loaded'
    #
    # id = np.load('identity_%s' % args.suffix)
    # print 'IDs loaded'
    #
    # dt = np.load('dt_%s' % args.suffix)
    # END temporary section

    no_comp = np.shape(pos)[1]
    nT = np.shape(pos)[2]

    dt = (trj_times[-1] - trj_times[0])/(nT - 1)
    f = open('dt_%s' % args.suffix, 'w')
    np.save(f, dt)
    f.close()

    f = open('%s' % args.trajfile, 'r')

    a = []
    for line in f:
        a.append(line)

    f.close()

    id = np.zeros([1, no_comp], dtype=object)

    count = 2

    for j in range(no_comp):
        atom = str.strip(a[count][10:15])
        id[0, j] = atom
        count += 1

    print 'id array made'
    f = open('identity_%s' % args.suffix, 'w')
    np.save(f, id)
    f.close()

    print 'IDs got and saved'


# Constants

kb = 1.38*10**-23  # Boltzmann constant [=] J/K
e = 1.602 * 10 ** -19  # elementary charge (C)
q = Atom_props.charge[args.ion] * e  # elementary charge (1.602e-19 C) * charge of ion
frontfrac = args.frontfrac
fracshow = args.fracshow
d = 3  # dimensionality


def nernst_einstein(D, D_std, C, C_std, T):
    NE_av = q**2*C*D/(kb*float(T))
    NE_error = NE_av*np.sqrt((D_std/D)**2 + (C_std/C)**2)
    return NE_av, NE_error

# if args.method == 'NE' and args.arrays == 'off' or args.method == 'B' and args.arrays == 'off':
# if args.arrays == 'off':
    # NOTE: we don't need to calculate concentration for CD but we do need z_1, z_2 and Lc which take the majority of
    # the time to calculate in the function below

    # print 'Calculating Concentration'
    # C_list = subprocess.check_output(["python", "Concentration.py", "-i", "%s" % args.trajfile, "-c", "%s" % args.ion,
    #                                   "-b", "%s" % args.buffer, "-l", "%s" % args.LC_type, "-s" "%s" % args.solv]).splitlines()
    # C_list = np.array([float(i) for i in C_list])
    # C, C_std, cross, Lc, z_2, z_1 = Concentration.conc('%s' % args.trajfile, '%s' % args.ion, '%s' % args.buffer,
    #                                                       '%s' % args.LC_type, '%s' % args.solv)
    # Lc = C_list[4] - C_list[3]
    # cross = C_list[2]
    # factor = 1*10**9  # number of nm in a m
    # conv = (Lc/factor)/(cross/(factor**2))
    # print 'Concentration Calculated'
    # Conc_props = np.array([C_list[0], C_list[1], Lc, conv, C_list[4], C_list[3]])

    # START Temporary section meant for debugging
    # Conc_props = np.load('C_props_%s' % args.suffix)
    # C = Conc_props[0]
    # C_std = Conc_props[1]
    # Lc = Conc_props[2]
    # conv = Conc_props[3]
    # z_2 = Conc_props[4]
    # z_1 = Conc_props[5]
    # print 'Properties loaded'
    # END temporary section meant for debugging


if args.method == 'NE' or args.method == 'B':

    if args.arrays == 'off':

        print 'Calculating Diffusivity'
        # Calculate diffusivity (m^2/s)
        _, _, _, D_av, D_std = Diffusivity.dconst(pos_ion, nT, args.nboot, frontfrac, fracshow, d, dt, args.nMC)
        D_tails = [D_av, D_std]  # Get it? D_tails, i.e. details about D

        f = open('D_tails_%s' % args.suffix, 'w')
        np.save(f, D_tails)
        f.close()

        print 'Diffusivity Saved'

    # D_tails = np.load('D_tails_array612ns')
    # D_av = D_tails[0]
    # D_std = D_tails[1]

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

plt.hist(ic_cds)
print "Ionic Conductivity: %s +/- %s S/m" % (ic_cd_avg, ic_cd_std)

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
