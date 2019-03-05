#!/usr/bin/env python

import matplotlib.pyplot as plt
import sqlite3 as sql
import numpy as np
import names

connection = sql.connect("../../../LLC_Membranes/timeseries/msd.db")
crsr = connection.cursor()
restrict_by_name = False 
MW = False  # plot mw of species
penalty = 0.25
wt = 5 

group = 'all'

restricted_groups = ['simple_alcohols', 'simple_plus_diols']
time_averaged = True 

if time_averaged:
	savename = '%s_%dwt_tamsds.pdf' % (group, wt)
else:
	savename = '%s_%dwt_emsds.pdf' % (group, wt)

query = "SELECT name, sim_length, MD_TAMSD, MD_TAMSD_CI_lower, MD_TAMSD_CI_upper, MD_MSD, MD_MSD_CI_lower, MD_MSD_CI_upper, mw from msd"

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

if wt == 10:
	query += " penalty = %s and" % penalty

query += " wt_water = %s" % wt
query += " and name != 'HOH'"
query += " ORDER BY MD_TAMSD DESC"

output = crsr.execute(query).fetchall()

labels = np.array([names.res_to_name[i[0]] for i in output], dtype=object)
sim_length = np.array([i[1] for i in output], dtype=float)
md_tamsd = np.array([i[2] for i in output])
md_tamsd_lower = -np.array([i[3] for i in output]) + md_tamsd
md_tamsd_upper = np.array([i[4] for i in output]) - md_tamsd
md_msd = np.array([i[5] for i in output], dtype=float)
md_msd_lower = -np.array([i[6] for i in output]) + md_msd
md_msd_upper = np.array([i[7] for i in output]) - md_msd
mw = np.array([i[8] for i in output])

# only need this block until all simulations are finished
md_msd *= (1000 / sim_length)
md_tamsd *= (400 / (sim_length * .4))

if time_averaged:

	ordered_md = np.argsort(md_tamsd)[::-1]
	md_tamsd = md_tamsd[ordered_md]
	md_tamsd_lower = md_tamsd_lower[ordered_md]
	md_tamsd_upper = md_tamsd_upper[ordered_md]
	mw = mw[ordered_md]

else:
	ordered_md = np.argsort(md_msd)[::-1]
	md_msd = md_msd[ordered_md]
	md_msd_lower = md_msd_lower[ordered_md]
	md_msd_upper = md_msd_upper[ordered_md]
	mw = mw[ordered_md]

labels = labels[ordered_md]
# end temporary block

fig, ax = plt.subplots(figsize=(6, 7))
bar_width = 0.8
opacity = 0.8
index = np.arange(len(labels), dtype=float)

if not MW:
	index += bar_width / 2

if time_averaged:
	ax.bar(index, md_tamsd, bar_width, alpha=opacity, yerr=(md_tamsd_lower, md_tamsd_upper), label='MD Simulated time-averaged MSD')
else:
	ax.bar(index, md_msd, bar_width, alpha=opacity, yerr=(md_msd_lower, md_msd_upper), label='MD Simulated MSD')

ax.set_ylabel('MSD ($nm^2$)', fontsize=14)
ax.tick_params(labelsize=14)
ax.set_xticks(index)# + bar_width/2)
ax.set_xticklabels(labels, fontsize=14)
plt.xticks(rotation=90)

if MW:
	ax1 = ax.twinx()
	ax1.bar(index + bar_width, mw, bar_width, alpha=opacity, label = 'Molecular Weight', color='red')
	ax1.tick_params(labelsize=14)
	ax1.set_ylabel('Molecular Weight (g/mol)', fontsize=14)

plt.xlim(-0.2, len(labels))
fig.tight_layout()
plt.savefig(savename)
plt.show()

