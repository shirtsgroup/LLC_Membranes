#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from LLC_Membranes.analysis.markov_state_dependent_dynamics import States
from LLC_Membranes.llclib import file_rw
from scipy.stats import levy_stable
import matplotlib.pyplot as plt
import numpy as np

res = 'URE'
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt/states.pl" % res
states = file_rw.load_object(path)

alpha, mu, sigma = states.fit_params[-1]
bins, edges = np.histogram(states.emissions[-1], bins=200, range=(-1.75, 1.75), density=True)
bin_width = edges[1] - edges[0]
centers = edges[:-1] + bin_width / 2
ratio = [bins[i] / levy_stable.pdf(x, alpha=alpha, beta=0, loc=mu, scale=sigma) for i, x in enumerate(centers)]

plt.figure(figsize=(10, 6))
plt.bar(centers, bins, bin_width, align='center', alpha=0.5, label='Empirical Disribution', color='xkcd:blue')
plt.plot(centers, levy_stable.pdf(centers, alpha=alpha, beta=0, loc=mu, scale=sigma), '--', lw=2,
         color='xkcd:orange', label='Analytical PDF')
plt.plot(centers, ratio, label='Ratio of Empirical:Analytical PDF', color='red', lw=2)
plt.plot(np.linspace(-1.75, 1.75, 3), np.ones([3]), '--', color='black', label='Exact agreement')
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=14)
plt.ylabel('Probability Density', fontsize=14)
plt.tight_layout()
plt.show()
