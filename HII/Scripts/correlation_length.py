#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import detect_peaks
from scipy.optimize import curve_fit


def exponential_decay(x, L):

    return 1 + np.exp(-x / L)


n = 20  # number of "layers"
spacing = 3.7  # distance between layers
L = 10  # correlation length
trials = 1000
nbins = 200
mpd = 10  # used for peak locating. Adjust if the wrong peaks are found. (higher for noisier data)

cov = np.zeros([n, n])  # initialize covariance matrix

z = np.linspace(0, (n-1)*spacing, n)
decay = np.exp(-z/L)  # decay of covariance

# decay[1:] += np.exp(-z[::-1][:-1]/L) # for periodicity

for i in range(z.shape[0]):
    cov[i, i:] = decay[:(n - i)]
    cov[i:, i] = decay[:(n - i)]

x = np.random.multivariate_normal(z, cov, trials)

# This method doesn't work because put all the data from each trial into a single array will just result in
# non-correlated gaussians

# x = x.flatten()

# duplicate periodically
# periodic = np.concatenate((x, x + n*spacing, x - n*spacing))

# d = np.zeros([x.shape[0], periodic.shape[0]])
# for i in range(x.shape[0]):
#    d[i, :] = np.abs(x[i] - periodic)

periodic = np.concatenate((x, x + n*spacing, x - n*spacing), axis=1)

d = np.zeros([trials, x.shape[1], periodic.shape[1]])  # pair-wise distances for each trial
h = np.zeros([trials, nbins])
for t in range(trials):
    for i in range(x.shape[1]):
        d[t, i, :] = np.abs(x[t, i] - periodic[t, :])
    h[t, :], edges = np.histogram(d[t, ...].flatten(), bins=nbins, range=[0.01, n*spacing/2])  # histogram trial

h = np.mean(h, axis=0)  # average histogram over all trials
centers = np.array([edges[i] + (edges[i + 1] - edges[i]) / 2 for i in range(h.shape[0])])  # centers of bins

h[0] = 0  # get rid of giant peak at zero
avg = np.mean(h)
h = np.array([h[i]/avg for i in range(len(h))]) # normalize

# find peaks to fit exponential to
peaks = detect_peaks.detect_peaks(h, mpd=mpd, show=False)  # adjust mpd if all peaks aren't found

if centers[peaks[0]] < spacing / 2:  # sometimes a peak is found where it shouldn't
    peaks = peaks[1:]

# plot averaged histogram
plt.bar(centers, h, centers[1] - centers[0])
# plot locations of peaks
plt.scatter(centers[peaks], h[peaks], marker='+', c='r', s=200, label='Peak locations')

# fit decaying exponential to peaks
p = np.array([L])  # initial guess at parameters
solp, cov_x = curve_fit(exponential_decay, centers[peaks], h[peaks], p)

# plot decaying exponential fit and decaying exponential with L that we were trying to match
plt.plot(centers, exponential_decay(centers, solp[0]), '--', c='black', label='Least squares fit')
plt.plot(centers, exponential_decay(centers, L), '--', c='blue', label='Theoretical')

print('Correlation length = %1.2f +/- %1.2f angstroms' % (solp[0], np.sqrt(cov_x[0, 0])))
plt.xlabel('Z distance separation (nm)', fontsize=14)
plt.ylabel('Count', fontsize=14)
plt.axes().tick_params(labelsize=14)
plt.legend(loc=0, prop={'size': 16})
plt.tight_layout()

# plt.figure()
# plt.imshow(cov)
plt.show()
