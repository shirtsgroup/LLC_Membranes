#!/usr/bin/env python

from LLC_Membranes.llclib import file_rw, fitting_functions, stats
from LLC_Membranes.analysis import Poly_fit
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

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
                vals = [maxes[i, c] for i, c in enumerate(choices)]
                b[boot], m[boot] = Poly_fit.poly_fit(np.log(L), np.log(vals), 1)[-1]

        return b, m	


def mfpt_analytical(passage_times, nbins, length, nboot, guess_params):

        passage_times /= 1000

        boot = np.zeros(nboot)

        for b in range(nboot):

            ptimes = np.random.choice(passage_times, size=passage_times.size, replace=True)

            hist, edges = np.histogram(ptimes, nbins, density=True)  # convert to microseconds
            bin_width = edges[1] - edges[0]
            bin_centers = np.array([i + bin_width / 2 for i in edges[:-1]])

            # very important to have a good guess. Might need to pass these unless I calculate MSD
            p0 = [length, guess_params[0] * length ** -guess_params[1], guess_params[2]]

            epsilon = 0.00001  # tolerance in parameter
            bounds = [(length - epsilon, 0, 0), (length + epsilon, np.inf, np.inf)]

            popt = curve_fit(fitting_functions.continuum_passage_time_distribution, bin_centers, hist, p0=p0,
                            bounds=bounds)[0]

            t = np.linspace(0.1, passage_times.max(), 1000)
            mfpt = quad(fitting_functions.continuum_ptime_distribution_expected_value, 0.1, np.inf,
                        args=(popt[0], popt[1], popt[2]))[0]

            boot[b] = mfpt
            #plt.plot(t, fitting_functions.continuum_passage_time_distribution(t, popt[0], popt[1], popt[2]), '--',
            #         color='black', lw=2)

            #plt.hist(ptimes, nbins, density=True)
            #plt.show()

        return boot

correlation = True 
root_dir = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11/"
nboot = 200 

#residues = ['URE', 'GCL', 'MET', 'ACH']
residues = np.array(['ACH', 'MET', 'URE', 'GCL'])
colors = {'URE':'xkcd:blue', 'GCL':'xkcd:orange', 'MET':'xkcd:green', 'ACH':'xkcd:magenta'}
names = {'URE': 'Urea', 'GCL': 'Ethylene Glycol', 'MET': 'Methanol', 'ACH': 'Acetic Acid'}
names2 = {'URE': 'Urea', 'GCL': 'Ethylene\nGlycol', 'MET': 'Methanol', 'ACH': 'Acetic\nAcid'}

guess_parameters = {'URE': [30, 1.5, 2.5], 'GCL': [30, 1.5, 2.5], 'MET': [np.exp(2.83), 1.7, 1], 'ACH': [np.exp(1.742), 1.6, 1]}
if correlation:
	Ls = {'URE': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'GCL': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'MET': [10, 15, 20, 25, 30, 35, 40, 45], 'ACH': [10, 15, 20, 25, 30, 35, 40, 45, 50]}
else:
	Ls = {'URE': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'GCL': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'MET': [10, 15, 20, 25, 30, 35, 40, 45, 50], 'ACH': [10, 15, 20, 25, 30, 35, 40, 45, 50]}

ndx = {10:0, 15:1, 20:2, 25:3, 30:4, 35:5, 40:6, 45:7, 50:8}
nbins = 100

########################
# Flux and MFPT curves #
########################

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

params = np.zeros([len(residues), nboot, 2]) 
labels = []
for j, r in enumerate(residues):
    print(r)
    if correlation:
        workdir = "%s/%s/10wt/mfpt/" % (root_dir, r)
    else:
        workdir = "%s/%s/10wt/brownian_mfpt" % (root_dir, r)

    L = Ls[r]

    mean = np.zeros(len(L))
    std = np.zeros(len(L))
    mfpts = np.zeros([len(L), nboot])

    for i, l in enumerate(L):
        print(l)
        passage_times = file_rw.load_object('%s/passage_times_%s_%d.pl' % (workdir, r, l))
        mfpts[i, :] = mfpt_analytical(passage_times, nbins, l, nboot, guess_parameters[r])
        mean[i], std[i] = mfpts[i, :].mean(), mfpts[i, :].std()

    fluxes = np.reciprocal(mfpts)
    flux_b, flux_m = bootstrap_m(nboot, fluxes, L)

    params[j, :, 0] = np.exp(flux_b)
    params[j, :, 1] = flux_m
    ax1.errorbar(L, fluxes.mean(axis=1), yerr=fluxes.std(axis=1), lw=2, elinewidth=2, capsize=5, capthick=2, color=colors[r])
    x = np.linspace(L[0], L[-1], 1000)
    if not correlation:
        ax1.plot(x, np.exp(flux_b.mean()) * x ** flux_m.mean() , '--', color=colors[r], label = r'%s ($\beta = %.2f \pm %.2f$)' % (names[r], flux_m.mean(), flux_m.std()), lw=2)
    else:
        ax1.plot(x, np.exp(flux_b.mean()) * x ** flux_m.mean() , '--', color=colors[r], label = '%s' % names[r], lw=2)

    b, m = bootstrap_m(nboot, mfpts, L)
    ax2.errorbar(L, mean, yerr=std, lw=2, elinewidth=2, capsize=5, capthick=2, color=colors[r], alpha=0.8)
    x = np.linspace(L[0], L[-1], 1000)
    #plt.plot(x, np.exp(b.mean()) * x ** m.mean() , '--', color=colors[r], label = r'%s ($\alpha = %.2f \pm %.2f, c = %.2f \pm %.2f $)' % (names[r], m.mean(), m.std(), np.exp(b.mean()), np.exp(b).std()), lw=2)
    #ax2.plot(x, np.exp(b.mean()) * x ** m.mean() , '--', color=colors[r], label = r'%s ($\beta = %.2f \pm %.2f$)' % (names[r], m.mean(), m.std()), lw=2)
    ax2.plot(x, np.exp(b.mean()) * x ** m.mean() , '--', color=colors[r], label = '%s' % names[r], lw=2)
    labels.append(Line2D([0], [0], color=colors[r], label='%s' % names[r], lw=2))

labels.append(Line2D([0], [0], color='black', linestyle='--', label=r'Fits to $cL^{-\beta}$', lw=2))

ax1.set_xlabel('Pore Length, L (nm)', fontsize=14)
ax1.set_ylabel('Flux ($\mu$s$^{-1}$)', fontsize=14)
ax1.tick_params(labelsize=14)
ax1.legend(handles=labels, fontsize=14)
fig1.tight_layout()

ax2.set_xlabel('Pore Length, L (nm)', fontsize=14)
ax2.set_ylabel('Mean First Passage Time ($\mu$s)', fontsize=14)
ax2.tick_params(labelsize=14)
ax2.legend(fontsize=14)
fig2.tight_layout()

if correlation:
    fig1.savefig('flux_curves.pdf')
    fig2.savefig('mfpt_curves.pdf')
else:
    fig1.savefig('flux_curves_brownian.pdf')
    fig2.savefig('mfpt_curves_brownian.pdf')
    plt.show()
    exit()

####################
# Selectivity Plot #
####################

#permutations = [[1, 2], [4, 2], [4, 1], [3, 1], [3, 2], [4, 3]]
permutations = [[4, 1], [4, 2], [3, 1], [3, 2], [4, 3], [1, 2]]
macroL = np.linspace(0.01, 100, 1000)

plt.figure()
for p in permutations:
        ratios = np.zeros([nboot, macroL.size])
        for b in range(nboot):
            choice1, choice2 =  np.random.choice(nboot, size=2)
            param1 = params[p[0] - 1, choice1, :]
            param2 = params[p[1] - 1, choice2, :]
            ratios[b, :] = (param1[0] * macroL ** param1[1]) / (param2[0] * macroL ** param2[1])
        error = stats.confidence_interval(ratios, 95)  # 95 % confidence interval
        mean = ratios.mean(axis=0)
        plt.plot(macroL, mean, label='%s:%s' % (names[residues[p[0] - 1]], names[residues[p[1] - 1]]), lw=2)
        plt.fill_between(macroL, mean + error[1, :], mean - error[0, :], alpha=0.5)

plt.legend(fontsize=12)
plt.xlabel('Pore length ($\mu m$)', fontsize=14)
plt.ylabel('Selectivity', fontsize=14)
plt.tick_params(labelsize=14)
plt.tight_layout()
#plt.savefig('selectivity.pdf')

#######################
# Parameter Bar Chart #
#######################

#fig, bar = plt.subplots()
#bar2 = bar.twinx()

bar_locations = np.array([1, 2, 3, 4])
bar_width1 = 0.8
bar_width2 = 0.4

c = np.array([p[:, 0].mean() for p in params])
c_err = np.array([p[:, 0].std() for p in params])

beta = np.array([p[:, 1].mean() for p in params])
beta_err = np.array([p[:, 1].std() for p in params])

#bar.bar(bar_locations - bar_width/2, c, bar_width, edgecolor='black', color='xkcd:blue', label='c', yerr=c_err)
#bar2.bar(bar_locations + bar_width/2, -beta, bar_width, edgecolor='black', color='xkcd:magenta', label=r'$\beta$', yerr=beta_err)

H = []
Herr = []
for r in residues:
    params = file_rw.load_object('%s/%s/10wt/forecast_%s_1state.pl' %(root_dir, r, r))
    H.append(np.mean(params.hurst_distribution))
    Herr.append(np.std(params.hurst_distribution))

H = np.array(H)
Herr = np.array(Herr)

#H = np.array([0.34, 0.30, 0.37, 0.40])

#for i, h in enumerate(H):
#    bar2.text(bar_locations[i] + bar_width/2, -beta[i] + 0.025, '%.2f' % h, fontsize=12, horizontalalignment='center', fontweight='bold', color='xkcd:green')

#hatch1 = mpatches.Patch(facecolor='xkcd:blue', label='c', edgecolor='black')
#hatch2 = mpatches.Patch(facecolor='xkcd:magenta', label=r'$\beta$', edgecolor='black')
#hatch3 = mpatches.Patch(facecolor='xkcd:green', label='H', edgecolor='black')

#bar.set_ylabel('c', fontsize=14)
#bar2.set_ylabel(r'$\beta$', fontsize=14)
#bar2.set_ylim(2, 3)
#plt.xticks(np.arange(1, 5), [names2[r] for r in residues])
#bar.tick_params(labelsize=14)
#bar2.tick_params(labelsize=14)
#plt.legend(handles=[hatch1, hatch2, hatch3], fontsize=14, loc='upper left')


fig1, bar1 = plt.subplots()

bar1.bar(bar_locations, c, bar_width1, edgecolor='black', color='xkcd:blue', label='c', yerr=c_err)
plt.xticks(np.arange(1, 5), [names2[r] for r in residues])
bar1.tick_params(labelsize=14)
bar1.set_ylabel('c', fontsize=14)

fig1.tight_layout()
fig1.savefig('c_parameters.pdf')

fig, ax = plt.subplots()
ax2 = ax.twinx()

reordered = [3, 2, 0, 1]

ax.bar(bar_locations - bar_width2/2, -beta[reordered], bar_width2, edgecolor='black', color='xkcd:blue', label=r'$\beta$', yerr=beta_err[reordered])
ax2.bar(bar_locations + bar_width2/2, H[reordered], bar_width2, edgecolor='black', color='xkcd:orange', label='$H$', yerr=Herr[reordered])

hatch1 = mpatches.Patch(facecolor='xkcd:blue', label=r'$\beta$', edgecolor='black')
hatch2 = mpatches.Patch(facecolor='xkcd:orange', label='$H$', edgecolor='black')
ax.legend(handles=[hatch1, hatch2], fontsize=14, loc='upper left', ncol=2)

ax.set_ylabel(r'$\beta$', fontsize=14)
ax.tick_params(labelsize=14)
ax2.set_ylabel('$H$', fontsize=14)
ax2.tick_params(labelsize=14)

ax.set_ylim(2.4, 2.8)
ax2.set_ylim(0.2, 0.45)
plt.xticks(np.arange(1, 5), [names2[r] for r in residues[reordered]])

fig.tight_layout()
fig.savefig('beta_parameters.pdf')

plt.show()
