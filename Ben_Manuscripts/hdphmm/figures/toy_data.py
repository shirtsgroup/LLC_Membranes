#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import file_rw
from hdphmm import generate_timeseries as gt
from hdphmm import timeseries as ts


def msd(traj, nboot=200, endshow=2000, dt=0.5, show=False, color='xkcd:blue', label=None):

    # Calculate MSD and plot
    msd = ts.msd(traj, 1)
    error = ts.bootstrap_msd(msd, nboot, confidence=68)[0]
    
    t = np.arange(endshow)*dt
    plt.plot(t, msd.mean(axis=1)[:endshow], lw=2, color=color, label=label)
    plt.fill_between(t, msd.mean(axis=1)[:endshow] + error[0, :endshow], msd.mean(axis=1)[:endshow] - error[1, :endshow], alpha=0.3, color=color)

    #for traj in range(msd.shape[1]):
    #    plt.plot(t, msd[:endshow, traj], color='black', lw=2)

    plt.tick_params(labelsize=14)
    plt.xlabel('Time (ns)', fontsize=14)
    plt.ylabel('Mean Squared Displacement (nm$^2$)', fontsize=14)

    plt.tight_layout()
    if show:
        plt.show()


def gentraj(params, ntraj, ndraws, alpha):
    
    trajectory_generator = gt.GenARData(params)
    trajectory_generator.gen_trajectory(ndraws, ntraj, bound_dimensions=[0], resample_T=False, alpha=alpha)

    return trajectory_generator.traj

ntraj = 24 
ndraws = 5000
alpha = 1000

#params = file_rw.load_object('saved_parameters/final_parameters_agglomerative_MET_diags_nsigma5_nA5_nr3.pl')
sigmas = np.array([[[0.10, -0.01], [-0.01, 0.11]],
                  [[0.01, 0], [0, 0.01]]])
As = np.array([[[0.25, -0.03], [-0.03, 0.25]],
              [[0.10, 0.01], [0.01, 0.10]]])

mu = np.zeros([24, 2])
mu[:, 0] = [0.7, 0.8, 0.9, 1.0, 0.6, 0.8, 0.9, 1.0, 1.1, 1.0, 0.7, 0.9, 1.9, 2.1, 2.3, 2.2, 2.4, 2.2, 2.5, 2.0, 2.1, 2.3, 1.8, 2.4]

Ts = [0.9, 0.98, 0.998]

sigma = np.zeros([24, 2, 2])
A = np.zeros([24, 1, 2, 2])
T = np.zeros([24, 24])

weights = [1, 5, 10]
for r in range(2):
    for a in range(2):
        for s in range(2):
            for t in range(3):
                ndx = r * 12 + a * 6 + s * 3 + t
                sigma[ndx, ...] = sigmas[s]
                A[ndx, 0, ...] = As[a]
                T[ndx, ndx] = Ts[t]

                #indices = np.random.choice([i for i in range(24) if i != ndx], size=10, replace=False)
                indices = np.array([i for i in range(24) if i != ndx])                
                total = 1 - Ts[t]
                weights = np.array([weights[i % 3] for i in indices])
                weights = weights / weights.sum()
                weights *= total

                T[ndx, indices] = weights


delete = [1, 2, 7, 8, 13, 14, 19, 20, 3, 9, 15, 21]
keep = [i for i in range(24) if i not in delete]

grid = tuple(np.meshgrid(keep, keep))
T = T[grid]

for i in range(len(keep)):
    ndx = [j for j in range(len(keep)) if j != i]
    tot = 1 - T[i, i]
    T[i, ndx] /= (T[i, ndx].sum() / tot)

print(T)
A = A[keep, ...]
sigma = sigma[keep, ...]
mu = mu[keep, :]

pi_init = np.ones(len(keep)) / len(keep)

params = {'T': T, 'A': A, 'sigma': sigma, 'pi_init': pi_init, 'mu': mu}

colors = ['xkcd:blue', 'xkcd:red', 'xkcd:green']

#for i, ntraj in enumerate([100, 200, 300]):
for i, alpha in enumerate([1000, 10000, 100000]):
    traj = gentraj(params, ntraj, ndraws, alpha)
    msd(traj, color=colors[i], label='alpha=%d' % alpha)

file_rw.save_object((traj, 0.5), 'toy_data.pl')

plt.legend(loc='upper left', fontsize=14)
#plt.savefig('msd_alpha_sensitivity.pdf')
plt.show()
exit()
msd(traj, show=True)

time = np.arange(ndraws) * 0.5
for t in range(ntraj):

    fig, ax = plt.subplots(2, 1, figsize=(12, 7))

    ax[0].plot(time, traj[:, t, 0])
    ax[1].plot(time, traj[:, t, 1])    

    plt.show()
