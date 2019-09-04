#!/usr/bin/env python

"""
This uses the Successive Random Additions method to generate approximations of fractional Brownian motion trajectory.
It kind of works, but does a bad job with autocorrelation.
"""

import numpy as np
import matplotlib.pyplot as plt
from LLC_Membranes.llclib import timeseries

tol = 1e-6
H = 0.3
n = 12
mean = 0
sigma_0 = 1
sigma_n = np.sqrt((1 - 2 ** (2*H - 2)) * (sigma_0 ** 2 / (2 ** (2*H))))
npoints = 2**n + 1

t = np.zeros([npoints])
t[0], t[-1] = np.random.normal(loc=mean, scale=sigma_0, size=2)
ndx = [0, int(2**n)]

for j in range(n):

    new_ndx = [int(i*2**(n - j - 1)) for i in range(2**(j + 1) + 1)]  # something like 2**(n+1) + 1

    interp = [i for i in new_ndx if i not in ndx]

    pts = t[ndx]
    for i in range(pts.size - 1):
        t[interp[i]] = np.mean([pts[i], pts[i + 1]])

    t[new_ndx] += np.random.normal(loc=mean, scale=sigma_n, size=2**(j + 1) + 1)

    sigma_n = np.sqrt(sigma_n ** 2 / (2 ** (2*H)))

    ndx = [int(i * 2 ** (n - j - 1)) for i in range(2**(j + 1) + 1)]

while sigma_n / sigma_0 > tol:

    t += np.random.normal(loc=mean, scale=sigma_n, size=t.size)

    sigma_n = np.sqrt(sigma_n ** 2 / (2 ** (2*H)))

t += np.random.normal(loc=mean, scale=sigma_n, size=t.size)
acf = timeseries.acf(t[1:] - t[:-1])

c = 0.5 * (2 ** (2*H) - 2)
print(c)
print(acf[1])
plt.plot(acf)
# plt.hist(t[1:] - t[:-1], bins=50)
plt.show()



