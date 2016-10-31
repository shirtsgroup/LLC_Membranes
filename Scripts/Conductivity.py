#!/usr/bin/python

"""
The purpose of this script is to calculate the ionic conductivity of a given LLC membrane using the Collective Diffusion
model described in the following paper:

    Y. Liu and F. Zhu "Collective diffusion model for ion conduction through microscopic channels," Biophysical Journal,
    vol 104, no. 2, pp. 368-376, Jan. 2013.

Ionic conductivity quantifies the movement of ions from site to site through defects in the crystal lattice of
ionic solids or aqueous solutions.

To calculate the ionic conductivity, we need to know the current of ions. For current through a channel we define the
current through the channel as:

    I = sum ( e_i * v_i / Lc )                                                                                (1)

where e and v are the charge and velocity of components i located in the channel. Lc is the channel length defined by
two planes at z1 and z2. Thus Lc = z2 - z1. Equation (1), on average, is the same as the current through the entire
system.

The charge transfer through the channel is defined as:

    Q(t) = integral[0, t](I(t)dt)                                                                             (2)

The equilibrium diffusion coefficient of Q, Dq, is given by:

    Dq = integral[0, infinity](dt*expected_value(I(0)*I(t)))                                                  (3)

The steady state current is given by:

    I_steady = (Dq * V) / (k_b * T)                                                                           (4)

Solving for I/V gives the system conductance:

    c = Dq / (k_b * T)                                                                                        (5)

This can be used to calculate the ionic conductance with an equilibrium simulation

Alternatively we can calculate conductance directly from current and voltage if an electric field is applied across the
membrane. To do so, calculate the average current using a time average of equation (1) and the potential difference
across the membrane (which can be found using gmx potential with a careful look at the output .xvg file)
"""

import argparse
from Get_Positions import get_positions
import Last_Frame
import Thickness
import Atom_props
import math
import matplotlib.pyplot as plt
import time
import Poly_fit
import numpy as np

start = time.time()
parser = argparse.ArgumentParser(description = 'Run Cylindricity script')

parser.add_argument('-f', '--trajfile', default='wiggle_traj.gro', help = 'Path to input file')
parser.add_argument('-i', '--ion', default='NA', help='Name of ion(s) being used to calculate ionic conductivity')
parser.add_argument('-l', '--LC_type', default='HII', help='Type of liquid crystal. Should match that defined in '
                                                           'LC_class.py')
parser.add_argument('-s', '--solv', default='no', help='Is the system solvate?')
parser.add_argument('-b', '--buffer', default='10', help='Distance into membrane from min and max where current '
                                                        'measurements will be made')
parser.add_argument('-T', '--temp', default=300, help='System Temperature, Kelvin')

args = parser.parse_args()

pos, vel, trj_times, box = get_positions('%s' % args.trajfile, '%s' % args.ion, '%s' % args.LC_type, '%s' % args.solv)

nT = len(trj_times)
Nbootstraps = 2000
frontfrac = 0
fracshow = 1
d = 1
dt = (trj_times[len(trj_times) - 1] - trj_times[0])/(len(trj_times) - 1)

no_ion = len(pos[0][0])
x = np.zeros([3, no_ion, len(trj_times)])
for i in range(3):
    for j in range(no_ion):
        for k in range(len(trj_times)):
            x[i, j, k] = pos[k][i][j]


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

MSD = np.zeros([nT], dtype=float)
MSDs = np.zeros([nT, no_ion], dtype=float)  # a set of MSDs per particle
for n in range(0, no_ion):
    MSDs[:, n] = msd_fft(x[:, n, :].T)
    MSD += MSDs[:, n]
MSD /= no_ion


eMSDs = np.zeros([nT, Nbootstraps], dtype=float)  # a set of MSDs per particle (time step?)
# Now, bootstrap over number of particles, assuming that the particles are sufficiently independent
for b in range(0, Nbootstraps):
    indices = np.random.randint(0, no_ion, no_ion)  # randomly select N of the particles with replacement
    # ^ makes a list of length Nparticles each with a random number from 0 to  Nparticles
    for n in range(0, no_ion):
        eMSDs[:, b] += MSDs[:, indices[n]]  # for this bootstrap trial b, add the MSDs of a randomly selected particle
        # to the particle at each time step
    eMSDs[:, b] /= no_ion  # Divide every timestep by Nparticles -- average the MSDs

# The result is Nbootstraps trajectories

limits = np.zeros([2, nT], dtype=float)  # a set of MSDs per particle
# now, let's determine a 95\% error bound for each tau (out of the
# Nbootstrapped MSD's, use that for the error bars
for t in range(0, nT):
    limits[0, t] = np.abs(np.percentile(eMSDs[t, :], 2.5) - MSD[t])
    limits[1, t] = np.abs(np.percentile(eMSDs[t, :], 97.5) - MSD[t])

# plot just the first half of the MSD (avoid the noisy part)
endMSD = int(np.floor(nT*fracshow))
startMSD = int(np.floor(nT*frontfrac))
print dt*startMSD
print dt*endMSD
# plot only 100 bars total
errorevery = int(np.ceil(fracshow*nT/100.0))
# Make an error bar plot
plt.errorbar(dt*np.array(range(0,endMSD)),MSD[:endMSD],yerr=[limits[0,:endMSD],limits[1,:endMSD]],errorevery=errorevery)
plt.ylabel('MSD')
plt.xlabel('time (ps)')
# I'm not sure why the error appears when plt.show is done.
y_fit2, r_squared2, std2, coeff2 = Poly_fit.poly_fit(dt*np.array(range(startMSD, endMSD)), MSD[startMSD:endMSD], 1)
print coeff2[1]/(2*d*100)  # slope of MSD plot equals 2 D. Divide by 0.01 to go from nm^2/ps to cm^2/s
#plt.plot(dt*np.array(range(startMSD, endMSD)), y_fit2)
plt.title('Stackoverflow/Shirts algorithm MSD results')
plt.show()

exit()
# UNITS
# pos [=] nm
# vel [=] nm/ps
# trj_times [=] ps
# box [=] nm

# Calculate z1 and z2 using z_max and z_min from Thickness.py

Last_Frame.extract_last_gro('%s' % args.trajfile)

# Calculate the thickness of the last frame

thick, z_max, z_min = Thickness.thickness('last_frame.gro')  # Name of output from Last_Frame.py

# Take a specified % off each z value to look at a smaller portion of the interior membrane to avoid end effects
# NOTE: in an infinite system without end effects, current can be measured along any Lc. I'm just being careful here

center = (z_max + z_min)/2
z_1 = center - (center - z_min)*(1 - float(args.buffer)/100.0)
z_2 = center + (z_max - center)*(1 - float(args.buffer)/100.0)

Lc = z_2 - z_1  # length of channel over which we will measure current

e = Atom_props.charge[args.ion] * 1.602 * 10 ** - 19  # charge of ions (1.602e-19 C)


def current(frame, positions, velocities, channel_length, charge, zmin, zmax):
    I = 0
    for i in range(0, len(positions[frame][0])):
        if zmin <= positions[frame][2][i] <= zmax:
            velocity = velocities[frame][2][i]
            I += velocity * charge / channel_length
    return I

currents = []
for i in range(0, len(trj_times)):
    currents.append(current(i, pos, vel, Lc, e, z_1, z_2))

time_step = trj_times[len(trj_times) - 1] / (len(trj_times) - 1)


def current_autocorrelation(I, steps):
    corr_sum = 0
    for i in range(0, steps + 1):
        corr_sum += I[0]*I[i]
    return corr_sum / (steps + 1)

def current_autocorrelation2(I, steps):
    itau = 0
    MSD = np.zeros([steps], dtype=float)  # find msd at each time step
    while itau < steps:
        for t in range(steps - itau):
            xo = I[t + itau] - I[t]
            MSD[itau] += xo**2
        MSD[itau] /= steps
        itau += 1
    return MSD

autocorrelation = []
for i in range(0, len(trj_times)):
    autocorrelation.append(current_autocorrelation(currents, i))

kb = 1.38*10**-23  # Boltzmann constant
# Dq = 1*10**12*sum(autocorrelation)/len(autocorrelation)  # [=] C^2 / s
# sigma = Dq/(kb*T)
# print sigma
# print Lc

# This is how the paper does it ...


def dQ(frame, positions, channel_length, zmax, zmin):
    q = 0
    for i in range(0, len(positions[frame][0])):  # find out how far each ion has moved in the z direction this frame
        if zmin <= positions[frame][2][i] <= zmax:  # check that the ion is in the channel over which we are measuring
            current_frame = positions[frame][2][i]  # If it, that is its position at the current frame
        elif positions[frame][2][i] > zmax:  # If the particle is above zmax, we are still potentially interested in it
            current_frame = zmax
        elif positions[frame][2][i] < zmin:  # If the particle is below zmin, we are still potentially interested in it
            current_frame = zmin
        if zmin <= positions[frame - 1][2][i] <= zmax:  # Do the same checks on the particle from the previous frame
            prev_frame = positions[frame - 1][2][i]
        elif positions[frame - 1][2][i] > zmax:
            prev_frame = zmax
        elif positions[frame - 1][2][i] < zmin:
            prev_frame = zmin
        if prev_frame == zmin and current_frame == zmin:
            pass
        if prev_frame == zmax and current_frame == zmax:
            pass
        displacement = current_frame - prev_frame
        q += displacement/channel_length
    return q

dq_frame = []
for i in range(1, len(trj_times)):
    dq_frame.append(dQ(i, pos, Lc, z_2, z_1))

q_sum = []
for i in range(0, len(dq_frame)):
    q_sum.append(sum(dq_frame[0:i+1]))
end = time.time()

print 'Total elapsed time: %s seconds' % (end - start)

y_fit, r_squared, std, coeff = Poly_fit.poly_fit(trj_times[1:len(trj_times)], q_sum, 1)

slope = coeff[1]  # the slope will be the second entry in the coeff vector output by poly_fit, i.e. the coefficient of x
Dq = slope / 2

# Cross Sectional Area of Box:
x_box_vect = math.sqrt(box[0][len(box[0]) - 1]**2 + box[3][len(box[3]) - 1]**2 + box[4][len(box[4]) - 1]**2)
y_box_vect = math.sqrt(box[1][len(box[1]) - 1]**2 + box[5][len(box[5]) - 1]**2 + box[6][len(box[6]) - 1]**2)
area = x_box_vect*y_box_vect

sigma = abs(1*10**12 * e**2 * Dq / (kb*float(args.temp))*(1*10**9*Lc/area))


print 'Ionic Conductivity: %s S/m' % sigma

plt.figure()
plt.title('Q(t), Ionic Conductivity = %.2e S/m' % sigma)
plt.plot(trj_times[1:len(trj_times)], q_sum, label='Raw Data')
plt.plot(trj_times[1:len(trj_times)], y_fit, label='Linear Fit')
plt.xlabel('Time (ps)')
plt.ylabel('Total Charge Transfer (C)')
# plt.text(1, 1, 'Ionic Conductivity: \n %s S' % sigma, ha='left', va='top')
plt.plot()
plt.figure()
test = current_autocorrelation2(currents, len(trj_times))
y_fit2, r_squared2, std2, coeff2 = Poly_fit.poly_fit(trj_times, test, 1)
plt.plot(trj_times, test)
plt.plot(trj_times, y_fit2)
plt.title('Test')
plt.show()


