#!/usr/bin/env python

from LLC_Membranes.analysis.coordination_number import System
from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt
import names
from LLC_Membranes.analysis import lifetime
# run: coordination_number.py -t PR_nojump.xtc -g berendsen.gro -r NA -rc res -at O / O3 O4

wt = 5 
path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11'
pickle = 'solute_NA_coordination.pl'
savename = 'all_NA_dwells.pdf'
ci = 0.95

residues = ['ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE']

mean_5wt = np.zeros([len(residues), 2])
ci_5wt = np.zeros([len(residues), 2])
mean_10wt = np.zeros([len(residues), 2])
ci_10wt = np.zeros([len(residues), 2])

try:
        loaded = np.load('na_lifetimes_ci_%.2f.npz' % ci)
        mean_5wt = loaded['mean_5wt']
        ci_5wt = loaded['ci_5wt']
        mean_10wt = loaded['mean_10wt']
        ci_10wt = loaded['ci_10wt']

except FileNotFoundError:
        for i, r in enumerate(residues):
                for wt in [10, 5]:
                        print(r, '%swt' % wt)
                        full_path = '%s/%s/%swt' % (path, r, wt)
                        life = lifetime.CoordinationLifetime('%s/PR_nojump.xtc' % full_path, '%s/berendsen.gro' % full_path, residue=r, coordinated_residue='NA', type='O')
                        life.calculate_lifetimes(ci=ci)
                        if wt == 5:
                                mean_5wt[i] = life.mean_lifetime
                                ci_5wt[i] = life.confidence
                        elif wt == 10:
                                mean_10wt[i] = life.mean_lifetime
                                ci_10wt[i] = life.confidence
                        print(life.mean_lifetime, life.confidence)

        np.savez_compressed('na_lifetimes_ci_%.2f' % ci, mean_5wt=mean_5wt, ci_5wt=ci_5wt, mean_10wt=mean_10wt, ci_10wt=ci_10wt)


labels = np.array([names.abbreviation[r] for r in residues])
colors = np.array([names.color_dict[r] for r in residues])

values = np.zeros([len(residues)])

#                       ACH    ACN    ATO    BUT    DMF    DMP    DMS    EAC    ETH    GCL    GLY    MET    PCB    PG     PR     RIB    SOH    TET    THF    URE
#values_5wt = np.array([3.808, 4.106, 5.935, 3.345, 5.143, 3.858, 4.080, 5.045, 3.372, 3.417, 3.792, 3.682, 3.134, 3.717, 3.476, 3.812, 3.471, 3.642, 3.310, 4.659]) / 2 # convert to ns

#                        ACH    ACN    ATO    BUT    DMF    DMP    DMS    EAC    ETH    GCL    GLY    MET    PCB    PG     PR     RIB    SOH    TET    THF    URE
#values_10wt = np.array([2.655, 3.425, 3.452, 1.989, 3.186, 2.050, 3.600, 3.981, 2.358, 2.337, 2.713, 2.009, 3.018, 2.760, 2.362, 3.122, 2.037, 2.735, 2.839, 3.398]) / 2

#ordered_10wt = ['MET', 'SOH', 'GCL', 'ACH', 'URE', 'ETH', 'ACN', 'TET', 'GLY', 'PR', 'DMS', 'PG', 'PCB', 'ATO', 'EAC', 'DMP', 'BUT', 'RIB', 'THF', 'DMF']
#ordered = [res.index(r) for r in ordered_10wt]
ordered = np.argsort(ci_10wt[:, 0])[::-1]
labels = labels[ordered]
colors = colors[ordered]

#fig, ax = plt.subplots(figsize=(6, 7))
fig, ax = plt.subplots()

#for i, r in enumerate(res):
#	print(r)
	
#	c = file_rw.load_object('%s/%s/%dwt/%s' % (path, r, wt, pickle))

#	n = [(c.distances[t] != 0).sum(1).mean() for t in range(len(c.distances))]
#	print(np.mean(n))

#	values[i] = 100 * np.mean(n)

index = np.arange(len(residues), dtype=float)
bar_width = 0.4
opacity =0.8
ax.bar(index, ci_10wt[ordered, 0], bar_width, label='10 wt % water', yerr=ci_10wt[ordered, 1])
ax.bar(index + bar_width, ci_5wt[ordered, 0], bar_width, alpha=opacity, label='5 wt % water', yerr=ci_5wt[ordered, 1])
ax.set_xticks(index + bar_width/2)
ax.set_xticklabels(labels, fontsize=14)
[x.set_color(colors[i]) for i, x in enumerate(plt.gca().get_xticklabels())]
plt.xticks(rotation=90)

plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.ylabel('95$^{th}$ Percentile Dwell Time (ns)', fontsize=14)
plt.gca().set_ylim(bottom=0.5)
plt.legend(fontsize=14)
fig.tight_layout()
plt.savefig(savename)
plt.show()

