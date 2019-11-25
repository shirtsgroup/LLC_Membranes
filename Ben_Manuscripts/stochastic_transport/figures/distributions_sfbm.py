#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from sympy import mpmath
from LLC_Membranes.analysis.sfbm_parameters import SFBMParameters
from LLC_Membranes.llclib import file_rw
from scipy.stats import levy_stable, norm
import matplotlib.pyplot as plt
import numpy as np

res = 'ACH'
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/%s/10wt/forecast_%s_1state.pl" % (res, res)
states = file_rw.load_object(path)

hops = ['Gaussian', 'Levy']
powerlaws = ['power', 'powerlaw_cutoff']

#hops = ['Gaussian']
#powerlaws = ['power']

powerlaw_savename = '%s_powerlaw.pdf' % res
#powerlaw_savename = '%s_powerlaw_%s' % (res, powerlaws[0])
hop_savename = 'gaussian_levy_comparison_anomalous_%s.pdf' % res
#hop_savename = '%s_hops_%s.pdf' % (hops[0], res)


# Correlation between hops
states.estimate_hurst(show=False, color='black', savename='%s_hop_acf.pdf' % res, plot_params=False)
plt.ylim(-0.5, 1)
plt.tight_layout()
plt.savefig('%s_hop_acf.pdf' % res)

fig, ax1 = plt.subplots()
ax2 = fig.add_axes([0.45, 0.3, 0.45, 0.4])
dwells = states.dwell_times[0]

nbins = 50
xmax = 200
xmin = 3
bins, edges = np.histogram(dwells, bins=nbins, range=(xmin, xmax), density=True)
bin_width = edges[1] - edges[0]
bin_centers = [i + bin_width / 2 for i in edges[:-1]]
ax1.bar(bin_centers, bins, bin_width, color='xkcd:blue', alpha=0.6)

ax2min = 200 
ax2max = 600
nbetween = 0
for i in dwells:
	if i < ax2max and i > ax2min:
		nbetween += 1

bins2, edges2 = np.histogram(dwells, bins=nbins, range=(ax2min, ax2max), density=True)

bin_width2 = edges2[1] - edges2[0]
bins2 *= (nbetween / len(dwells))   # get the density right
bin_centers2 = [i + bin_width2 / 2 for i in edges2[:-1]]
ax2.bar(bin_centers2, bins2, bin_width2, color='xkcd:blue', alpha=0.6)

x = np.linspace(xmin, xmax, 500)
x2 = np.linspace(ax2min, ax2max, 500)
for dwell in powerlaws:
	states.fit_distributions(plot=False, show=False, dwell_distribution=dwell, nboot=10)	
	
	if dwell == 'powerlaw_cutoff':
		mean_alpha = np.mean([p[0] for p in states.dwell_parameters[0]]) + 1
		mean_lambda = np.mean([p[1] for p in states.dwell_parameters[0]])
		y = x ** -mean_alpha * np.exp(-mean_lambda * x)	
		y2 = x2 ** -mean_alpha * np.exp(-mean_lambda * x2)	
		C = mean_lambda ** (1 - mean_alpha) / float(mpmath.gammainc(1 - mean_alpha, mean_lambda * xmin))
		y *= C
		y2 *= C
		xstart = np.argmin(np.abs(y - bins[0]))
		ax1.plot(x[xstart:], y[xstart:], '--', color='xkcd:orange', label='Truncated Power Law', lw=2)
		ax2.plot(x2, y2, '--', color='xkcd:orange', lw=2)
	elif dwell == 'power':
		mean_alpha = np.mean(states.dwell_parameters[0])
		C = (mean_alpha) * xmin **(mean_alpha)
		y = x ** -(1 + mean_alpha)
		y2 = x2 ** -(1 + mean_alpha)
		y *= C
		y2 *= C
		xstart = np.argmin(np.abs(y - bins[0]))
		ax1.plot(x[xstart:], y[xstart:], '--', color='xkcd:green', label='Power Law', lw=2)
		ax2.plot(x2, y2, '--', color='xkcd:green', label='Power Law', lw=2)

ax1.legend(fontsize=14)
ax1.set_xlabel('Dwell time (ns)', fontsize=14)
ax1.set_ylabel('Probability', fontsize=14)
ax1.tick_params(labelsize=14)
ax2.tick_params(labelsize=14)
#ax2.set_yticks([]) # turn of y-axis
plt.tight_layout()
plt.savefig(powerlaw_savename)

all_hops = []

plt.figure()
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
plt.xlabel('Hop Length (nm)', fontsize=14)
plt.xlim(-1.5, 1.5)
plt.tight_layout()
plt.savefig(hop_savename)
plt.show()
