#!/usr/bin/env python

from LLC_Membranes.analysis.coordination_number import System
from LLC_Membranes.analysis import coordination_number
from LLC_Membranes.llclib import file_rw
import numpy as np
import matplotlib.pyplot as plt
import names

# run: coordination_number.py -t PR_nojump.xtc -g berendsen.gro -r NA -rc res -at O / O3 O4

head_groups = True
path = '/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11'
pickle = 'solute_NA_coordination.pl'
savename = 'all_solutes_NA_coordination.pdf'

residues = ['GLY', 'ACH', 'ACN', 'ATO', 'BUT', 'DMF', 'DMP', 'DMS', 'EAC', 'ETH', 'GCL', 'GLY', 'MET', 'PCB', 'PG', 'PR', 'RIB', 'SOH', 'TET', 'THF', 'URE']

labels = np.array([names.abbreviation[r] for r in residues])
colors = np.array([names.color_dict[r] for r in residues])

nsolute = 24

try:
	loaded = np.load('na_coord.npz')
	n_5wt = loaded['n_5wt']
	n_10wt = loaded['n_10wt']

except FileNotFoundError:

	n_5wt = np.zeros([len(residues), 2])
	n_10wt = np.zeros([len(residues), 2])

	for wt in [10, 5]:
		for i, r in enumerate(residues):
			loc = '%s/%s/%swt' % (path, r, wt)
			print(loc)
			# unnecessary since looking at unique solute-sodum associations (i.e. a solute coordinated to two sodium ions doesn't count boi)
			#sys = file_rw.load_object('%s/NA_O_coordination.pl' % loc)
			#if r in ['ACH', 'EAC', 'PCB']:
			#	if r == 'ACH':
			#		atoms = ['O1']
			#	elif r == 'EAC':
			#		atoms = ['O1']
			#	elif r == 'PCB':
			#		atoms = ['O2']
			#	sys = coordination_number.System('%s/PR_nojump.xtc' % loc, '%s/berendsen.gro' % loc, residue=r, coordinated_residue='NA', atoms=atoms)
			#else:
			sys = coordination_number.System('%s/PR_nojump.xtc' % loc, '%s/berendsen.gro' % loc, residue=r, coordinated_residue='NA', type='O')

			sys.distance_search(cut=0.25)
			sys.n_coordinated(plot=False)

			# narrow down to unique sodium - oxygen interactions
			nT = sys.t.n_frames
			ncoord = np.zeros([nT, nsolute])
			n_O = sys.ncoord.shape[1] // nsolute  # number of oxygens per residue 
			
			for t in range(nT):
				for o in range(nsolute):
					if sys.ncoord[t, n_O*o:n_O*(o + 1)].sum() > 0:
						ncoord[t, o] = 1
			
			# bootstrap -- assumes independent solutes
			nboot = 200
			n = []
			for b in range(nboot):
				choices = np.random.randint(nsolute, size=nsolute)
				nperframe = ncoord[:, choices].sum(axis=1)
				n.append(100 * np.mean(nperframe)/ nsolute)

			if wt == 5:
				n_5wt[i] = [np.mean(n), np.std(n)]
			if wt == 10:
				n_10wt[i] = [np.mean(n), np.std(n)]

	np.savez_compressed('na_coord.npz', n_5wt=n_5wt, n_10wt=n_10wt)		

#                         ACH    ACN     ATO    BUT     DMF     DMP    DMS    EAC    ETH    GCL      GLY     MET      PCB    PG        PR       RIB        SOH      TET       THF        URE
#values_10wt = np.array([0.2900, 0.477, 0.264, 0.0915, 0.2904, 0.0694, 0.350, 0.2700, 0.122, 0.1615, 0.1546, 0.09302, 0.3000, 0.1697, 0.124979, 0.131129, 0.104375, 0.13825, 0.1469375, 0.60833]) * 100

#                         ACH    ACN     ATO    BUT        DMF       DMP     DMS          EAC     ETH         GCL      GLY        MET      PCB    PG        PR       RIB        SOH      TET       THF        URE
#values_5wt = np.array([0.5500, 0.8336, 0.5258, 0.223125, 0.6885625, 0.23908, 0.7241541, 0.380, 0.270125, 0.303604, 0.3466944, 0.2658333, 0.51, 0.33555, 0.32333, 0.260141, 0.27404, 0.25372, 0.25725, 0.782625]) * 100

#ordered_10wt = ['MET', 'SOH', 'GCL', 'ACH', 'URE', 'ETH', 'ACN', 'TET', 'GLY', 'PR', 'DMS', 'PG', 'PCB', 'ATO', 'EAC', 'DMP', 'BUT', 'RIB', 'THF', 'DMF']
#ordered = [res.index(r) for r in ordered_10wt]

ordered = np.argsort(values_10wt)[::-1]
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

index = np.arange(len(res), dtype=float)
bar_width = 0.4
opacity =0.8
ax.bar(index, values_5wt[ordered], bar_width, alpha=opacity, label='5 wt % water')
ax.bar(index + bar_width, values_10wt[ordered], bar_width, label='10 wt % water')
ax.set_xticks(index + bar_width/2)
ax.set_xticklabels(labels, fontsize=14)
[x.set_color(colors[i]) for i, x in enumerate(plt.gca().get_xticklabels())]
plt.xticks(rotation=90)

plt.gcf().get_axes()[0].tick_params(labelsize=14)
plt.ylabel('Percent of Solute Coordinated to Sodium', fontsize=14)
plt.legend(fontsize=14)
fig.tight_layout()
plt.savefig(savename)
plt.show()

