#!/usr/bin/env python

from LLC_Membranes.llclib import file_rw
from hdphmm.generate_timeseries import GenARData
import numpy as np
import matplotlib.pyplot as plt

#np.random.seed(10)

params = file_rw.load_object('/home/ben/github/hdphmm/notebooks/saved_parameters/final_parameters_agglomerative_MET_eigs_combined_30.pl')

#print(params['T'])

trajectory_generator = GenARData(params=params)
trajectory_generator.gen_trajectory(4808, 24, bound_dimensions=[0])

hops = np.array(trajectory_generator.hops)[:, 1].tolist()
dwells = trajectory_generator.dwells

fig, ax = plt.subplots(1, 2, figsize=(12, 5))

# dwells
bins, edges = np.histogram(dwells, 25, range=(0, 150), density=True)
#bins3, edges3 = np.histogram(pre_cluster_dwells, 25, range=(0, 150), density=True)

bin_width = edges[1] - edges[0]
centers = [i + bin_width/2 for i in edges[:-1]]

ax[0].plot(centers, bins, lw=2, label='MD Data')
#ax[0].plot(centers, bins3, lw=2, label='IHMM Simulation Data')
# hops

bins_hop, edges_hop = np.histogram(hops, 51, range=(-2, 2), density=True)
#bins3_hop, edges3_hop = np.histogram(pre_cluster_hops, 51, range=(-2, 2), density=True)

bin_width_hop = edges_hop[1] - edges_hop[0]
centers_hop = [i + bin_width_hop/2 for i in edges_hop[:-1]]

#std = np.std(pre_cluster_hops)  # fit hops to Gaussian
#alpha, _, scale = levy.fit_levy(pre_cluster_hops, beta=0)[0].x  # fit hops to levy stable

#x = np.linspace(-2, 2, 500)
#levy = stats.levy_stable.pdf(x, alpha, 0, loc=0, scale=scale)
#normal = stats.norm.pdf(x, scale=std)

# ax[1].plot(x, normal, lw=2, label='$\mathcal{N}$(0, %.2f)' %std)
# ax[1].plot(x, levy, lw=2, label=r'$\mathcal{L}$(%.2f, %.2f)' %(alpha, scale * np.sqrt(2)))
#print('alpha = %.2f' % alpha)
ax[1].plot(centers_hop, bins_hop, lw=2, label='MD Data')
#ax[1].plot(centers_hop, bins3_hop, lw=2, label='IHMM Simulation Data')

ax[0].set_xlabel('Dwell Time (steps)', fontsize=14)
ax[0].set_ylabel('Probability', fontsize=14)
ax[0].tick_params(labelsize=14)
ax[0].legend(fontsize=14)

# file_rw.save_object({'hop_centers': centers_hop, 'dwell_centers': centers, 'x': x, 'levy': levy, 'normal': normal,
#                      'raw_hops': pre_cluster_hops, 'raw_dwells': pre_cluster_dwells, 'MD_dwells': dwells,
#                      'MD_hops': hops}, 'dwell_hops_ar%d.pl' % ihmm[0].order)

ax[1].set_xlabel('Jump length (nm)', fontsize=14)
ax[1].set_ylabel('Probability', fontsize=14)
ax[1].tick_params(labelsize=14)
ax[1].legend(fontsize=14)

fig.tight_layout()
plt.show()

#bins_hop, edges_hop = np.histogram(hops, 51, range=(-2, 2), density=True)

#bin_width_hop = edges_hop[1] - edges_hop[0]
#centers_hop = [i + bin_width_hop/2 for i in edges_hop[:-1]]

#plt.plot(centers_hop, bins_hop, lw=2, label='MD Data')

#plt.plot(trajectory_generator.traj[:, 0, 1])
#plt.show()
