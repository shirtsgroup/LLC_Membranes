#!/usr/bin/env python

import matplotlib.pyplot as plt
import sqlite3 as sql
import numpy as np
import names

connection = sql.connect("../../../timeseries/msd.db")
crsr = connection.cursor()
restrict_by_name = False 
penalty = 0.5

group = 'all'

restricted_groups = ['simple_alcohols', 'simple_plus_diols']
time_averaged = True 

if time_averaged:
	savename = '%s_tamsds.pdf' % group
else:
	savename = '%s_emsds.pdf' % group

query = "SELECT name, mw from msd"

if restrict_by_name:
	#restricted_names = ['MET', 'ETH', 'PR', 'BUT']
	#restricted_names = ['ETH', 'MET', 'BUT', 'PR', 'GCL', 'GLY', 'PG']
	#restricted_names = ['ETH', 'PR', 'BUT', 'RIB', 'TET']
	#restricted_names = ['GCL', 'SOH', 'GLY', 'BUT', 'DMP']
	restricted_names = ['ATO', 'ACH', 'DMS', 'URE', 'ACN']
	#restricted_names = ['EAC', 'THF', 'PCB']
	query += ' WHERE ('
	for i, n in enumerate(restricted_names):
		query += " name = '%s'" % n
		if i < len(restricted_names) - 1:
			query += ' or'
	query += ') and'
else:
	query += ' WHERE'

query += " penalty = %s" % penalty
query += " ORDER BY MD_TAMSD DESC"

output = crsr.execute(query).fetchall()

labels = np.array([names.res_to_name[i[0]] for i in output], dtype=object)
mw = np.array([i[1] for i in output])

fig, ax = plt.subplots(figsize=(12, 7))
bar_width = 0.8
opacity = 0.8
index = np.arange(len(labels))

ax.bar(index + bar_width/2, mw, bar_width, alpha=opacity, label='MD Simulated time-averaged MSD')

ax.set_ylabel('MSD ($nm^2$)', fontsize=14)
ax.tick_params(labelsize=14)
ax.set_xticks(index + bar_width/2)
ax.set_xticklabels(labels, fontsize=14)
plt.xticks(rotation=90)
plt.xlim(-0.2, len(labels))
fig.tight_layout()
plt.savefig(savename)
plt.show()

