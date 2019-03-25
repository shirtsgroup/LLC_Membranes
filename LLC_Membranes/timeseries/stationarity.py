#!/usr/bin/env python

"""
Testing for stationarity.
Based on Chapter 3 of Time Series Analysis by James Hamilton

A process is "weakly stationary" or "covariance-stationary" if:
E(Y_t) = mu     for all t
E(Y_t - mu)(Y_(t-j) - mu) = gamma_j = 0 for j != 0 and sigma^2 for j = 0

j : number of observations
mu : mean
gamma_j : covariance
sigma^2 : variance
"""

import numpy as np
import matplotlib.pyplot as plt
from statsmodels.tsa.stattools import adfuller


def autocovariance(joint_distribution):

    """
    E(yt)(ytlag)
    In words: expected value of yt times ytlag. They are not independent so you can't assume it equals E(yt)*E(ytlag)
    """

    observations = joint_distribution.shape[1]
    autocov = np.zeros([observations])

    for j in range(observations):

        yt_expected_value = joint_distribution[:, -1] - joint_distribution[:, -1].mean()
        ytlag_expected_value = joint_distribution[:, -j] - joint_distribution[:, -j].mean()

        autocov[j] = (yt_expected_value * ytlag_expected_value).mean()

    return autocov


# Generate data
nrealizations = 1000  # number of times to repeat sample generation (number of vectors in x p.45)
observations = 10000  # number of observations in each realization (y_t, y_(t-1), y_(t-2) ... y_(t-observations) p.45)
standard_deviation = 0.5  # standard deviation of normal distribution
mean = 0  # average of standard deviation

joint_distribution = np.zeros([nrealizations, observations])

for i in range(nrealizations):
    joint_distribution[i, :] = np.random.normal(loc=mean, scale=standard_deviation, size=observations)

autocov = autocovariance(joint_distribution)

means = np.zeros([observations])

for j in range(observations):
    means[j] = np.mean(joint_distribution[:, j])

fig, ax = plt.subplots(1, 2)
t = np.linspace(0, observations, observations)
ax[0].plot(t, autocov)
ax[0].set_xlabel('time (ps)')
ax[0].set_title('$\gamma_j$')
ax[1].plot(t, means)
ax[1].set_title('$\mu$')
ax[1].set_xlabel('time (ps)')

result = adfuller(joint_distribution[0, :])  # want p-value as low as possible if stationary

print(joint_distribution[0, :])
print('ADF Statistic: %f' % result[0])
print('p-value: %f' % result[1])
print('Critical Values:')
for key, value in result[4].items():
    print('\t%s: %.3f' % (key, value))

# plt.show()