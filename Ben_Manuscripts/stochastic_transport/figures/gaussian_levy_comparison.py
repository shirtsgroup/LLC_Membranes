#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from LLC_Membranes.analysis.markov_state_dependent_dynamics import States
from LLC_Membranes.llclib import file_rw
from scipy.stats import levy_stable, norm
import matplotlib.pyplot as plt
import numpy as np

res = 'URE'
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt/states.pl" % res
states = file_rw.load_object(path)
 
# Levy stable fit
alpha, mu, sigma = states.fit_params[-1]
mu_gaus = np.mean(states.emissions[-1])
sigma_gaus = np.std(states.emissions[-1])

bins, edges = np.histogram(states.emissions[-1], bins=200, range=(-1.75, 1.75), density=True)
bin_width = edges[1] - edges[0]
centers = edges[:-1] + bin_width / 2
ratio = [bins[i] / levy_stable.pdf(x, alpha=alpha, beta=0, loc=mu, scale=sigma) for i, x in enumerate(centers)]

plt.figure(figsize=(8, 5))
plt.bar(centers, bins, bin_width, align='center', alpha=0.5, label='Empirical Disribution', color='xkcd:blue')

x = np.linspace(-0.75, 0.75, 500)
plt.plot(x, levy_stable.pdf(x, alpha=alpha, beta=0, loc=mu, scale=sigma), '--', lw=2,
         color='xkcd:orange', label=r'L$\'{e}$vy Stable PDF')
plt.plot(x, norm.pdf(x, loc=mu_gaus, scale=sigma_gaus), '--', lw=2, color='xkcd:green', label='Gaussian PDF')
plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=14)
plt.ylabel('Probability Density', fontsize=14)
plt.xlabel('Hop length (nm)', fontsize=14)
plt.xlim(-0.75, 0.75)
plt.tight_layout()
plt.savefig('gaussian_levy_comparison.pdf')
plt.show()
