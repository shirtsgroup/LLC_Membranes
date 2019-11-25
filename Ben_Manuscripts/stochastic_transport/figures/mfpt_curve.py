#!/usr/bin/env python

from LLC_Membranes.llclib import file_rw
from LLC_Membranes.analysis import Poly_fit
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def growth(x, c, alpha):

        return c * x ** alpha

def bootstrap(nboot, passage_times, nbins=100, percentile=99.5):
	
        max_pt = np.zeros(nboot)
        limits = (passage_times.min(), np.percentile(passage_times, percentile))

        bar_heights = np.zeros(nbins)

        for b in range(nboot):

            ptimes = np.random.choice(passage_times, size=passage_times.size, replace=True)

            n, edges = np.histogram(ptimes, bins=nbins, range=limits, density=True)
            bin_width = edges[1] - edges[0]
            bin_centers = [x + bin_width / 2 for x in edges[:-1]]

            max_pt[b] = bin_centers[np.argmax(n)]
            bar_heights += n

        bar_heights /= nboot

        return max_pt

def bootstrap_m(nboot, maxes, L):

        b = np.zeros([nboot])
        m = np.zeros([nboot])
        for boot in range(nboot):
                choices = np.random.choice(maxes.shape[1], size=maxes.shape[0], replace=True)
                vals = [maxes[i, c] / 1000 for i, c in enumerate(choices)]
                b[boot], m[boot] = Poly_fit.poly_fit(np.log(L), np.log(vals), 1)[-1]

        return b, m	

correlation = False  

root_dir = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/"

#residues = ['URE', 'GCL', 'MET', 'ACH']
residues = ['ACH', 'MET', 'URE', 'GCL']
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
names = {'URE': 'Urea', 'GCL': 'Ethylene Glycol', 'MET': 'Methanol', 'ACH': 'Acetic Acid'}

if correlation:
	Ls = {'URE': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'GCL': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'MET': [10, 15, 20, 25, 30, 35, 40, 45], 'ACH': [10, 15, 20, 25, 30, 35, 40, 45]}
else:
	Ls = {'URE': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'GCL': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'MET': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'ACH': [10, 15, 20, 25, 30, 35, 40, 45, 50]}

ndx = {10:0, 15:1, 20:2, 25:3, 30:4, 35:5, 40:6, 45:7, 50:8}

for r in residues:

    if correlation:
        workdir = "%s/%s/10wt/mfpt/" % (root_dir, r)
    else:
        workdir = "%s/%s/10wt/brownian_mfpt" % (root_dir, r)

    nboot = 200
    L = Ls[r]

    mean = np.zeros(len(L))
    std = np.zeros(len(L))
    maxes = np.zeros([len(L), nboot])

    for i, l in enumerate(L):
        passage_times = file_rw.load_object('%s/passage_times_%s_%d.pl' % (workdir, r, l))
        maxes[i, :] = bootstrap(nboot, passage_times, nbins=100)
        mean[i], std[i] = maxes[i, :].mean(), maxes[i, :].std()

    b, m = bootstrap_m(nboot, maxes, L)
    print(r, b.mean(), m.mean())

    plt.errorbar(L, mean / 1000, yerr=std / 1000, lw=2, elinewidth=2, capsize=5, capthick=2, color=colors[r], alpha=0.8)
    x = np.linspace(L[0], L[-1], 1000)
    #plt.plot(x, np.exp(b.mean()) * x ** m.mean() , '--', color=colors[r], label = r'%s ($\alpha = %.2f \pm %.2f, c = %.2f \pm %.2f $)' % (names[r], m.mean(), m.std(), np.exp(b.mean()), np.exp(b).std()), lw=2)
    plt.plot(x, np.exp(b.mean()) * x ** m.mean() , '--', color=colors[r], label = r'%s ($\alpha = %.2f \pm %.2f$)' % (names[r], m.mean(), m.std()), lw=2)

plt.xlabel('Pore Length, L (nm)', fontsize=14)
plt.ylabel('Mean First Passage Time ($\mu$s)', fontsize=14)
plt.tick_params(labelsize=14)
plt.legend(fontsize=12)
plt.tight_layout()

if correlation:
    plt.savefig('mfpt_curves.pdf')
else:
    plt.savefig('mfpt_curves_brownian.pdf')

plt.show()
