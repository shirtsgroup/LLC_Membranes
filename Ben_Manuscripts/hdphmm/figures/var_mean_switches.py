#!/usr/bin/env python

import numpy as np
from hdphmm.generate_timeseries import GenARData
from hdphmm.utils import timeseries as ts
import matplotlib.pyplot as plt

T = np.array([[0.995, 0.005], [0.005, 0.995]])

phis = np.array([[[0.5, 0.1], [0.2, -0.3]],
    [[0.5, 0.1], [0.2, -0.3]]])[:, np.newaxis, ...]

cov = np.array([[[1, 0.005], [0.005, 1]],
                [[0.01, 0.005], [0.005, 0.02]]])

mu = np.array([[0, 5],
               [0, 1]])

pi_init = np.array([0.5, 0.5])

bound = [1]

params = {'A': phis, 'sigma':cov, 'mu':mu, 'T':T, 'pi_init':pi_init}

data = GenARData(params=params, dim=2, transition_matrix=T, phis=phis, order=1, cov=cov)

data.gen_trajectory(10000, 1, bound_dimensions=[1])

fig, ax = plt.subplots(2, 1)
ax[0].plot(data.traj[:, 0, 0])
ax[1].plot(data.traj[:, 0, 1])
plt.show()

exit()

# The following is implemented in generate_timeseries.py
ntraj = 1
order = 1
ndraws = 10000
traj = np.zeros([ndraws + order, ntraj, 2])

state_labels = [0, 1]
initial_state = np.random.choice(state_labels, p=pi_init)

# draw state sequence
state_sequence = np.zeros([ndraws], dtype=int)
state_sequence[0] = initial_state
for d in range(1, ndraws):
    previous_state = state_sequence[d - 1]
    state_sequence[d] = np.random.choice(state_labels, p=T[previous_state, :])

switch_points = ts.switch_points(state_sequence)
states = state_sequence[switch_points[:-1]]

tot_steps = 0

for i, sp in enumerate(switch_points[:-1]):

    nsteps = switch_points[i + 1] - switch_points[i]
    state = states[i]

    subtraj = np.zeros([nsteps + order, 2])
    for d in range(order, nsteps + order):
        subtraj[d, :] = sum([phis[state, 0, ...] @ subtraj[d - (i + 1), :] for i in range(order)])
        subtraj[d, :] += np.random.multivariate_normal(np.zeros(2), cov[state, ...])

    shift = subtraj[order, :] - traj[tot_steps, 0, :]
    shift[bound] = subtraj[order:, bound].mean() - mu[state, bound]
    traj[(tot_steps + order):(tot_steps + order + nsteps), 0, :] = subtraj[order:, :] - shift 

    tot_steps += nsteps

#state = 0
#nsteps = 20
#subtraj = np.zeros([nsteps + order, 2])
#for d in range(order, nsteps + order):
#    subtraj[d, :] = sum([phis[state, 0, ...] @ subtraj[d - (i + 1), :] for i in range(order)])
#    subtraj[d, :] += np.random.multivariate_normal(mu[state, :], cov[state, ...])

#shift = subtraj[order, :] - traj[tot_steps, 0, :]
#traj[(tot_steps + order):(tot_steps + order + nsteps), 0, :] = subtraj[order:, :] - shift 

#tot_steps += nsteps
#
#state = 1
#nsteps = 50
#subtraj = np.zeros([nsteps + order, 2])
#for d in range(order, nsteps + order):
#    subtraj[d, :] = sum([phis[state, 0, ...] @ subtraj[d - (i + 1), :] for i in range(order)])
#    subtraj[d, :] += np.random.multivariate_normal(mu[state, :], cov[state, ...])

#shift = subtraj[order, :] - traj[tot_steps, 0, :]
#traj[(tot_steps + order):(tot_steps + order + nsteps), 0, :] = subtraj[order:, :] - shift 

#tot_steps += nsteps

fig, ax = plt.subplots(2, 1)
#ax[0].plot(np.arange(100), traj[order:(100 + order), 0, 0])
#ax[1].plot(np.arange(100), traj[order:(100 + order), 0, 1])
#ax[0].plot(100 + np.arange(20), subtraj[order:(20 + order), 0])
#ax[1].plot(100 + np.arange(20), subtraj[order:(20 + order), 1])
#ax[0].plot(100 + np.arange(20), subtraj[order:(20 + order), 0] - shift[0])
#ax[1].plot(100 + np.arange(20), subtraj[order:(20 + order), 1] - shift[1])
ax[0].plot(traj[order:tot_steps + order, 0, 0])
ax[1].plot(traj[order:tot_steps + order, 0, 1])
#ax[0].plot(state_sequence)
#ax[1].plot(state_sequence)

plt.show()
exit()

