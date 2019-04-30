#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import names
from LLC_Membranes.analysis import hbonds
from LLC_Membranes.analysis.hbonds import System
from LLC_Membranes.llclib import file_rw

residues = ['ACH', 'URE', 'ACN', 'ATO']
ndx = [0, 13, 1]
names = [names.res_to_name[i] for i in residues]
nsolute = 24
nboot = 200

loaded = np.load('unique_hbonds.npz')
hg = loaded['n_10wt'][ndx, :]
hg = np.concatenate((hg, np.array([[0, 0]])))

ndx = [0, 19, 1, 2]
loaded = np.load('na_coord.npz')
na = loaded['n_10wt'][ndx, :]

#hg = (100 / 24) * np.array([12.7, 1.84, 1.12, 0])
#don = (100 / 24) * np.array([6.6, 2.6, 1.2, 0])
#acc = (100 / 24) * np.array([1.45, 2.2, 1.6, 0.5])
# na = (100 / 24) * np.array([.28, .52, .45, .26]) * 24

path = "/home/bcoscia/Documents/Gromacs/Transport/NaGA3C11"
acc = np.zeros([len(residues), nboot])
don = np.zeros([len(residues), nboot])

for i, r in enumerate(residues):

        loc = '%s/%s/10wt' % (path, r)
        print(loc)
        acceptors = file_rw.load_object('%s/%s_HOH_acc_hbonds.pl' %(loc, r))  # hbonds.pl should have been run looking for hbonds between head group and monomer only (i.e. no water)
        
        res_no, nres = acceptors.number_residues(r)

        acceptor_matrix = np.zeros([acceptors.t.n_frames, nsolute], dtype=bool)

        for frame in range(acceptors.t.n_frames):
                res_numbers = []
                for j in acceptors.hbonds[frame][2]:
                    try:
                        res_numbers.append(res_no[j])
                    except KeyError:
                        pass
#               res_numbers = [res_no[i] for i in acceptors.hbonds[frame][0]]
                acceptor_matrix[frame, np.unique(res_numbers).astype(int)] = True

        for b in range(nboot):
                ndx = np.random.randint(0, nsolute, size=nsolute)
                acc[i, b] = acceptor_matrix[:, ndx].sum(axis=1).mean()

        donors = file_rw.load_object('%s/%s_HOH_donors_hbonds.pl' %(loc, r))
        res_no, nres = donors.number_residues(r)
        
        donor_matrix = np.zeros([donors.t.n_frames, nsolute], dtype=bool)
        for frame in range(donors.t.n_frames):
                res_numbers = [res_no[int(i)] for i in donors.hbonds[frame][0]]
                donor_matrix[frame, np.unique(res_numbers).astype(int)] = True

        for b in range(nboot):
                ndx = np.random.randint(0, nsolute, size=nsolute)
                don[i, b] = donor_matrix[:, ndx].sum(axis=1).mean()

don /= nsolute  # normalize
don *= 100  # convert to percent
acc /= nsolute
acc *= 100

index = np.arange(len(names))

fig, ax = plt.subplots()

bar_width = 0.2

ax.set_ylabel('% Solutes Involved Per Frame', fontsize=14)
ax.set_xlabel('   ')
ax.tick_params(labelsize=14)
ax.set_xticks(index + 2*bar_width)
ax.set_xticklabels(names, fontsize=14)

opacity = 0.7

ax.bar(index, hg[:, 0], bar_width, alpha=opacity, label='hbonds w/ head groups', yerr=hg[:, 1])
ax.bar(index + bar_width, don.mean(axis=1), bar_width, alpha=opacity, label='hbonds donated to water', yerr=don.std(axis=1))
ax.bar(index + 2*bar_width, acc.mean(axis=1).mean(), bar_width, alpha=opacity, label='hbonds accepted from water', yerr=acc.std(axis=1))
ax.bar(index + 3*bar_width, na[:, 0], bar_width, alpha=opacity, label='coordinated to sodium', yerr=na[:, 1])

plt.ylim(0, 85)
plt.legend(fontsize=14)
plt.tight_layout()
plt.savefig('ketone_hbonds.pdf')
plt.show()
