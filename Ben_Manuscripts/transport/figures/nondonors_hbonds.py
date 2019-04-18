#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import names
from LLC_Membranes.analysis import hbonds

#names = [names.res_to_name[i] for i in ['THF', 'PCB', 'EAC', 'DMF']]
names = ["Tetrahydrofuran", "Propylene\nCarbonate", "Ethyl\nAcetate", "Dimethyl\nFormamide"]
residues = ['THF', 'PCB', 'EAC', 'DMF']
#hg = [12.7, 1.84, 1.12, 0]
#donors = [6.6, 2.6, 1.2, 0]
acc = 100 * np.array([0.4475, 1.5, .65, .765]) / 24

loaded = np.load('na_coord.npz')

ndx = [18, 12, 7, 4]  # indices of residues based on residues list in all_coord.py
na = loaded['n_10wt'][ndx]
#na = 100 * np.array([.157, .292, .27, .28]) #* 24

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
wt = 10
nsolute = 24
acc = np.zeros([len(residues), 2])

try: 

	loaded = np.load('acceptors.npz')
	acc = loaded['acc']

except FileNotFoundError:

	for i, r in enumerate(residues):

		loc = '%s/%s/%swt' % (path, r, wt)

		sys = hbonds.System('%s/PR_nojump.xtc' % loc, '%s/berendsen.gro' % loc)
		sys.set_eligible(r, 'all')
		sys.set_eligible('HOH', 'all', donors_only=True)	
		sys.identify_hbonds(0.35, 30)

		#sys = file_rw.load_object('%s/hbonds.pl' %loc)
	
		n_unique_acceptors = []
		acceptors = np.zeros([sys.t.n_frames, nsolute], dtype=bool)

		res_no, nres = sys.number_residues(r)

		for frame in range(sys.t.n_frames):
			res_numbers = [res_no[i] for i in sys.hbonds[frame][2]]
			n_unique_acceptors.append(len(np.unique(res_numbers)))
			acceptors[frame, np.unique(res_numbers).astype(int)] = True

		# bootstrap
		n_unique_donors = np.array(n_unique_acceptors)
		nboot = 200
		n = []
		for b in range(nboot):
			choice = np.random.randint(nsolute, size=nsolute)
			nperframe = acceptors[:, choice].sum(axis=1)
			n.append(100*np.mean(nperframe) / nsolute)

		acc[i] = [np.mean(n), np.std(n)]

		np.savez_compressed('acceptors.npz', acc=acc)


index = np.arange(len(names))

fig, ax = plt.subplots()

bar_width = 0.4

ax.set_ylabel('Percent Solutes Involved Per Frame', fontsize=14)
ax.set_xlabel('   ')
ax.tick_params(labelsize=14)
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(names, fontsize=14)

opacity = 0.7

ax.bar(index, na[:, 0], bar_width, alpha=opacity, label='Coordinated to sodium', yerr=na[:, 1])
#ax.bar(index + bar_width, donors, bar_width, alpha=opacity, label='hbonds donated to water')
ax.bar(index + bar_width, acc[:, 0], bar_width, alpha=opacity, label='Hbonds accepted from water', yerr=acc[:, 1])
#ax.bar(index + 3*bar_width, na, bar_width, alpha=opacity, label='coordinated to sodium')

plt.ylim(0, 40)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('nondonor_hbonds.pdf')
plt.show()
