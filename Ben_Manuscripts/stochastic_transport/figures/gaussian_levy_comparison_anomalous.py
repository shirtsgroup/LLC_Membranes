#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
from LLC_Membranes.llclib import file_rw
from scipy.stats import levy_stable, norm
import matplotlib.pyplot as plt
import numpy as np

res = 'GCL'
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt/forecast_%s_1state.pl" % (res, res)
states = file_rw.load_object(path)

hops = ['Gaussian', 'Levy']

all_hops = []

for h in states.hop_lengths:
        all_hops += h[0]

plt.hist(all_hops, bins=50, color='xkcd:blue', alpha=0.6, density=True)
x = np.linspace(-1.5, 1.5, 500)

for h in hops:
	states.fit_distributions(plot=False, show=False, hop_distribution=h, nboot=10)

	if h == 'Gaussian':
                mean_sigma = np.mean([p[1] for p in states.hop_parameters[0]])
                mean_mu = np.mean([p[0] for p in states.hop_parameters[0]])
                plt.plot(x, norm.pdf(x, loc=mean_mu, scale=mean_sigma), '--', lw=2, color='xkcd:orange', label='Gaussian PDF')
	elif h == 'Levy':
                mean_alpha = np.mean([p[0] for p in states.hop_parameters[0]])
                mean_sigma = np.mean([p[2] for p in states.hop_parameters[0]])
                mean_mu = np.mean([p[1] for p in states.hop_parameters[0]])
                plt.plot(x, levy_stable.pdf(x, alpha=mean_alpha, beta=0, loc=mean_mu, scale=mean_sigma), '--', color='xkcd:green', label=r'L$\'{e}$vy PDF', lw=2)

plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.legend(fontsize=14)
plt.ylabel('Probability Density', fontsize=14)
plt.xlim(-1.5, 1.5)
plt.tight_layout()
plt.savefig('gaussian_levy_comparison_anomalous_%s.pdf' % res)
plt.show()
