#!/usr/bin/env python

import numpy as np
from LLC_Membranes.llclib import file_rw
import matplotlib.pyplot as plt
import levy
from scipy import stats

orders = [1, 2, 3, 4, 10]

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

for i, r in enumerate(orders):

    data = file_rw.load_object('/home/ben/github/hdphmm/notebooks/dwell_hops_ar%d.pl' % r)

    dwells = data['MD_dwells']
    pre_cluster_dwells = data['raw_dwells']
    hops = data['MD_hops']
    pre_cluster_hops = np.array(data['raw_hops'])

    pre_cluster_hops = pre_cluster_hops[np.where(pre_cluster_hops < 2)[0]]
    pre_cluster_hops = pre_cluster_hops[np.where(pre_cluster_hops > -2)[0]]

    bins, edges = np.histogram(dwells, 25, range=(0, 150), density=True)
    bins3, edges3 = np.histogram(pre_cluster_dwells, 25, range=(0, 150), density=True)

    bin_width = edges[1] - edges[0]
    centers = [i + bin_width/2 for i in edges[:-1]]

    if i == 0:
        ax[0].plot(centers, bins, lw=2, label='MD Data', color='black')

    ax[0].plot(centers, bins3, lw=2, label='AR(%d)' % r)
    # hops

    bins_hop, edges_hop = np.histogram(hops, 50, range=(-2, 2), density=True)
    bins3_hop, edges3_hop = np.histogram(pre_cluster_hops, 50, range=(-2, 2), density=True)

    bin_width_hop = edges_hop[1] - edges_hop[0]
    centers_hop = [i + bin_width_hop/2 for i in edges_hop[:-1]]

    std = np.std(pre_cluster_hops)  # fit hops to Gaussian
    alpha, _, scale = levy.fit_levy(pre_cluster_hops, beta=0)[0].x  # fit hops to levy stable

    x = np.linspace(-2, 2, 500)
    levy_ = stats.levy_stable.pdf(x, alpha, 0, loc=0, scale=scale)
    normal = stats.norm.pdf(x, scale=std)

    # ax[1].plot(x, normal, lw=2, label='$\mathcal{N}$(0, %.2f)' %std)
    # ax[1].plot(x, levy_, lw=2, label=r'$\mathcal{L}$(%.2f, %.2f)' %(alpha, scale * np.sqrt(2)))
    print('alpha = %.2f' % alpha)

    if i == 0:
        ax[1].plot(centers_hop, bins_hop, lw=2, label='MD Data', color='black')
    ax[1].plot(centers_hop, bins3_hop, lw=2, label='AR(%d)' % r)

ax[0].set_xlabel('Dwell Time (steps)', fontsize=14)
ax[0].set_ylabel('Probability', fontsize=14)
ax[0].tick_params(labelsize=14)
ax[0].legend(fontsize=14)

ax[1].set_xlabel('Jump length (nm)', fontsize=14)
ax[1].set_ylabel('Probability', fontsize=14)
ax[1].tick_params(labelsize=14)

ax[1].legend(fontsize=14)

fig.tight_layout()
plt.show()
