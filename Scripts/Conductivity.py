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
import Diffusivity

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
parser.add_argument('-B', '--nboot', default=2000, help='Number of bootstrap trials to be run')

args = parser.parse_args()

pos, vel, trj_times, box = get_positions('%s' % args.trajfile, '%s' % args.ion, '%s' % args.LC_type, '%s' % args.solv)

##### UNITS #####
# pos [=] nm
# vel [=] nm/ps
# trj_times [=] ps
# box [=] nm

no_ion = len(pos[0, :, 0])
nT = len(trj_times)
Nbootstraps = 2000
frontfrac = 0
fracshow = .3
d = 1
dt = (trj_times[len(trj_times) - 1] - trj_times[0])/(len(trj_times) - 1)

# Calculate diffusivity (cm^2/ps)
D = Diffusivity.dconst(pos, nT, Nbootstraps, frontfrac, fracshow, d, dt)[0]

exit()


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


