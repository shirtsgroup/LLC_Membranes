#!/usr/bin/env python
"""
Rescaled implementation follows:
https://en.wikipedia.org/wiki/Hurst_exponent

Clarified further here: http://bearcave.com/misl/misl_tech/wavelets/hurst/index.html
"""

import numpy as np
import fbm
import matplotlib.pyplot as plt
import tqdm

H = 0.4  # hurst parameter used to generate data
N = int(2 ** 15)  # number of points in time series
min_block = 8

Xt = fbm.FBM(N, H).fbm()[1:]
block_size = N

R = []
S = []
ns = []

count = 0
while block_size > min_block:

    ns.append(block_size)

    nsegments = int(N / block_size)

    for i in range(nsegments):

        partial_series = Xt[i*block_size:(i + 1)*block_size]
        m = partial_series.mean()
        Yt = partial_series - m  # mean-adjusted series

        try:
            R[count] += max(partial_series) - min(partial_series)
            S[count] += np.std(partial_series)
        except IndexError:
            R.append(max(partial_series) - min(partial_series))
            S.append(np.std(partial_series))

    R[count] /= nsegments
    S[count] /= nsegments

    block_size = int(block_size / 2)

    count += 1

plt.plot(np.log(ns), np.log(np.array(R) / np.array(S)))
plt.show()




