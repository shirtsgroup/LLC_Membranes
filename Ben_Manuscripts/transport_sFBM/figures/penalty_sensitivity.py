#!/usr/bin/env python

import sqlite3 as sql
import names
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(12, 7))
opacity = 0.8

connection = sql.connect("../../../timeseries/msd.db")
crsr = connection.cursor()

penalties = [0.25, 0.5, 0.62, 0.75, 0.88, 1.0]
colors = ['blue', 'red', 'green', 'xkcd:orange', 'purple', 'pink']

bar_width = 1 / (len(penalties) + 1)


query = "SELECT name, sigma, alpha, hurst, python_MSD, python_MSD_CI_lower, python_MSD_CI_upper from msd WHERE penalty = %.2f ORDER BY python_MSD DESC"  % penalties[0]
output = crsr.execute(query).fetchall()

labels = [names.res_to_name[i[0]] for i in output]
sigma = np.array([i[1] for i in output])
alpha = np.array([i[2] for i in output])
hurst = np.array([i[3] for i in output])
python_msd = np.array([i[4] for i in output])
python_msd_lower = -np.array([i[5] for i in output]) + python_msd
python_msd_upper = np.array([i[6] for i in output]) - python_msd
index = np.arange(len(labels))

ax.bar(index, python_msd, bar_width, alpha=opacity, color=colors[0], yerr=(python_msd_lower, python_msd_upper), label='Penalty = %.2f' % penalties[0])


for i, penalty in enumerate(penalties[1:]):

	query = "SELECT name, sigma, alpha, hurst, python_MSD, python_MSD_CI_lower, python_MSD_CI_upper from msd WHERE penalty = %.2f ORDER BY python_MSD DESC"  % penalty
	output = crsr.execute(query).fetchall()

	new_labels = [names.res_to_name[i[0]] for i in output]#, dtype=object)
	indices = np.array([new_labels.index(i) for i in labels])

	#new_labels = np.array(new_labels, dtype=object)

	#for i in range(20):
	#	print(labels[i], new_labels[indices][i])
	#exit()

	sigma = np.array([i[1] for i in output])
	alpha = np.array([i[2] for i in output])
	hurst = np.array([i[3] for i in output])
	python_msd = np.array([i[4] for i in output])
	python_msd_lower = -np.array([i[5] for i in output]) + python_msd
	python_msd_upper = np.array([i[6] for i in output])- python_msd

	ax.bar(index + (i+1)*bar_width, python_msd[indices], bar_width, alpha=opacity, color=colors[i+1], yerr=(python_msd_lower[indices], python_msd_upper[indices]), label='Penalty = %.2f' % penalty)

ax.set_ylabel('MSD ($nm^2$)', fontsize=14)
ax.tick_params(labelsize=14)
ax.set_xticks(index + bar_width/2)
ax.set_xticklabels(labels, fontsize=14)
ax.legend(fontsize=14)
plt.xticks(rotation=90)
fig.tight_layout()
plt.show()

