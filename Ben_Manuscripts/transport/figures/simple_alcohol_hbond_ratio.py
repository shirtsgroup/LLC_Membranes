#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import names
from LLC_Membranes.llclib import file_rw
from LLC_Membranes.analysis.hbonds import System
# generating using varia
# hbonds.py -t PR_nojump.xtc -g berendsen.gro -r HII BUT -a O O1 O2 O3 O4 -a all
# O O1 O2 are ethers. O3 O4 are carboxylate oxygens.

ndx = [7, 4, 9, 2]  # indices of alcohol entries in unique_hbonds array (see unique_hbonds.py)
frames = 2000
nsolute = 24
nboot = 200

loaded = np.load('unique_hbonds.npz')
total_hbonds = loaded['n_10wt'][ndx, :]
print(total_hbonds[:, 0], total_hbonds[:, 1])
resnames = ['MET', 'ETH', 'PR', 'BUT']

# incomplete trajectories and stricter geometric criteria
#ratio = [2.3, 3.1, 3.47, 5.51]
#total_hbonds = [8819, 9762, 8857, 8651]
#ratio = [1.21781242618, 1.62801302932, 1.98014411529, 2.65940967134]
#total_hbonds = 100 * np.array([18776, 20170, 18661, 17481]) / (frames * 24)
#total_hbonds = np.array([43.2, 42.2, 39.5, 38.9])

ether_oxygens = ["O", "O1", "O2"]
hg_oxygens = ["O3", "O4"]


ratio = np.zeros([len(resnames), nboot])
path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
for i, r in enumerate(resnames):

	loc = '%s/%s/10wt' % (path, r)
	print(loc)
	sys = file_rw.load_object('%s/hbonds.pl' %loc)  # hbonds.pl should have been run looking for hbonds between head group and monomer only (i.e. no water)

	# Find when an hbond is donated to an ether oxygen and when one is donated to a head group
	ethers = []
	hg = []
	for a in sys.t.topology.atoms:
		if a.name in ether_oxygens and a.residue.name == "HII":
			ethers.append(a.index)
		elif a.name in hg_oxygens and a.residue.name == "HII":
			hg.append(a.index)

	# need these arrays for bootstrapping
	ether_hb = np.zeros([sys.t.n_frames, nsolute], dtype=bool)
	hg_hb = np.zeros_like(ether_hb)

	res_no, nres = sys.number_residues(r)

	for frame in range(sys.t.n_frames):
		res_numbers = [res_no[i] for i in sys.hbonds[frame][0]]
		hbs = sys.hbonds[frame][2]  # hydrogen bond acceptors; sys.hbonds: [D, H, A, angle]
		hb_eth = [res_numbers[i] for i in range(len(hbs)) if hbs[i] in ethers]
		hb_hg = [res_numbers[i] for i in range(len(hbs)) if hbs[i] in hg]
		ether_hb[frame, hb_eth] = True
		hg_hb[frame, hb_hg] = True

	for b in range(nboot):
		n = np.random.randint(0, nsolute, size=nsolute)
		ratio[i, b] = hg_hb[:, n].sum() / ether_hb[:, n].sum()

index = np.arange(len(resnames))
opacity=0.75
bar_width = 0.4

fig, ax = plt.subplots()
 
ax.bar(index, ratio.mean(axis=1), bar_width, color='xkcd:blue', alpha=opacity, yerr=ratio.std(axis=1))
#ax.bar(index, ratio, bar_width, color='xkcd:blue', alpha=opacity)
ax.tick_params(axis='both', labelsize=14)
ax.set_ylabel('Carboxylate : Ether Hydrogen Bonds', fontsize=14, color='xkcd:blue')
ax.tick_params(axis='y', labelcolor='xkcd:blue')
ax2 = ax.twinx()
ax2.bar(index + bar_width, total_hbonds[:, 0], bar_width, color='xkcd:red', alpha=opacity, yerr=total_hbonds[:, 1])
ax2.tick_params(axis='both', labelsize=14, labelcolor='xkcd:red')
ax2.set_ylabel('% Solutes Hydrogen Bonded Per Frame', fontsize=14, color='xkcd:red')

ax.set_xticks(index + bar_width/2)
ax.set_xticklabels([names.res_to_name[r] for r in resnames], fontsize=14)
plt.tight_layout()
plt.savefig('simple_alcohol_hbonds.pdf')
plt.show()
