#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

alpha = 1000
x = np.zeros(10)
x[0] = 0.99
x[1:] = (1 - x[0]) / 9

colors = ['xkcd:blue', 'xkcd:red', 'xkcd:green']

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

for i, alpha in enumerate([500, 1000, 10000]):

    draws = np.random.dirichlet(x * alpha, size=100000)[:, 0]


    ax[0].hist(draws, 50, color=colors[i], label='alpha=%d' % alpha, density=True)
    ax[1].hist(1 / (1 - draws), 50, color=colors[i], label='alpha=%d' % alpha, density=True)

ax[0].legend(fontsize=14)
ax[1].legend(fontsize=14)

ax[0].set_xlabel('Self Transition Probability', fontsize=14)
ax[0].set_ylabel('Probability', fontsize=14)

ax[1].set_xlabel('Expected Dwell Time', fontsize=14)

ax[0].tick_params(labelsize=14)
ax[1].tick_params(labelsize=14)

ax[1].set_xlim(0, 750)

plt.tight_layout()
plt.show()
