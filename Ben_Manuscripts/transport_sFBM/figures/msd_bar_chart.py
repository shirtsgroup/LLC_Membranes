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

query = "SELECT name, python_MSD, python_MSD_CI_lower, python_MSD_CI_upper, MD_MSD, MD_MSD_CI_lower, MD_MSD_CI_upper, sigma, alpha, hurst, sim_length, MD_TAMSD, MD_TAMSD_CI_lower, MD_TAMSD_CI_upper, python_TAMSD, python_TAMSD_CI_lower, python_TAMSD_CI_upper from msd"

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
query += " ORDER BY MD_MSD DESC"

output = crsr.execute(query).fetchall()

labels = np.array([names.res_to_name[i[0]] for i in output], dtype=object)
python_msd = np.array([i[1] for i in output])
python_msd_lower = -np.array([i[2] for i in output]) + python_msd
python_msd_upper = np.array([i[3] for i in output]) - python_msd
md_msd = np.array([i[4] for i in output])
md_msd_lower = -np.array([i[5] for i in output]) + md_msd
md_msd_upper = np.array([i[6] for i in output]) - md_msd
sigma = np.array([i[7] for i in output])
alpha = np.array([i[8] for i in output])
hurst = np.array([i[9] for i in output])
sim_length = np.array([i[10] for i in output])
md_tamsd = np.array([i[11] for i in output])
md_tamsd_lower = -np.array([i[12] for i in output]) + md_tamsd
md_tamsd_upper = np.array([i[13] for i in output]) - md_tamsd
python_tamsd = np.array([i[14] for i in output])
python_tamsd_lower = -np.array([i[15] for i in output]) + python_tamsd
python_tamsd_upper = np.array([i[16] for i in output]) - python_tamsd

# only need this block until all simulations are finished
md_msd *= (1000 / sim_length)
md_tamsd *= (400 / (sim_length * .4))

if time_averaged:

	ordered_md = np.argsort(md_tamsd)[::-1]
	python_tamsd = python_tamsd[ordered_md]
	python_tamsd_lower = python_tamsd_lower[ordered_md]
	python_tamsd_upper = python_tamsd_upper[ordered_md]
	md_tamsd = md_tamsd[ordered_md]
	md_tamsd_lower = md_tamsd_lower[ordered_md]
	md_tamsd_upper = md_tamsd_upper[ordered_md]

else:
	ordered_md = np.argsort(md_msd)[::-1]
	python_msd = python_msd[ordered_md]
	python_msd_lower = python_msd_lower[ordered_md]
	python_msd_upper = python_msd_upper[ordered_md]
	md_msd = md_msd[ordered_md]
	md_msd_lower = md_msd_lower[ordered_md]
	md_msd_upper = md_msd_upper[ordered_md]

labels = labels[ordered_md]
# end temporary block

fig, ax = plt.subplots(figsize=(12, 7))
bar_width = 0.4
opacity = 0.8
index = np.arange(len(labels))

if time_averaged:
	rect2 = ax.bar(index , md_tamsd, bar_width, alpha=opacity, color='red', yerr=(md_tamsd_lower, md_tamsd_upper), label='MD Simulated time-averaged MSD')
	rects1 = ax.bar(index + bar_width, python_tamsd, bar_width, alpha=opacity, color='blue', yerr=(python_tamsd_lower, python_tamsd_upper), label='sFBM simulated Time-averaged MSD')
else:
	rect2 = ax.bar(index, md_msd, bar_width, alpha=opacity, color='red', yerr=(md_msd_lower, md_msd_upper), label='MD Simulated MSD')
	rects1 = ax.bar(index + bar_width, python_msd, bar_width, alpha=opacity, color='blue', yerr=(python_msd_lower, python_msd_upper), label='sFBM Simulated MSD')

#ax2 = ax.twinx()
#ax2.bar(index + bar_width, sigma, bar_width, alpha=opacity, color='green', label='Hurst Parameter')
#ax2.set_ylabel('Hurst parameter', fontsize=14)
#ax2.tick_params(labelsize=14)

#ax.set_xlabel('Molecule')
ax.set_ylabel('MSD ($nm^2$)', fontsize=14)
ax.tick_params(labelsize=14)
ax.set_xticks(index + bar_width/2)
ax.set_xticklabels(labels, fontsize=14)
ax.legend(fontsize=14)
plt.xticks(rotation=90)
fig.tight_layout()
plt.savefig(savename)
plt.show()

