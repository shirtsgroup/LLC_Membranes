#!/usr/bin/env python

import sqlite3 as sql
import names
import matplotlib.pyplot as plt
import numpy as np

time_averaged = True 

if time_averaged:
	mdmsd = "MD_TAMSD"
	pymsd = "python_TAMSD"
else:
	mdmsd = "MD_MSD"
	pymsd = "python_MSD"

fig, ax = plt.subplots(figsize=(12, 7))
opacity = 0.8

connection = sql.connect("../../../timeseries/msd.db")
crsr = connection.cursor()

penalties = [0.25, 0.5, 0.62, 0.75, 0.88, 1.0]
colors = ['blue', 'red', 'green', 'xkcd:orange', 'purple', 'pink']
names = ['URE', 'GCL', 'ACH', 'ETH']
bar_width = (1 / (len(penalties) + 1)) 

MD_query = "SELECT DISTINCT %s from msd WHERE" % mdmsd
for i, n in enumerate(names):
        MD_query += " name = '%s' and %s is not NULL" % (n, mdmsd)
        if i < len(names) - 1:
                MD_query += ' or'

MD_query += " ORDER BY name"
output = crsr.execute(MD_query).fetchall()
MD_MSDs = np.array([i[0] for i in output])

sFBM_MSDs = np.zeros([len(names), len(penalties)])

for j, pen in enumerate(penalties):

	query = "SELECT %s from msd WHERE (penalty = %.2f) AND ("  % (pymsd, pen)

	for i, n in enumerate(names):
        	query += " name = '%s'" % n
	        if i < len(names) - 1:
        	        query += ' or'

	query += " ) ORDER BY name"

	output = crsr.execute(query).fetchall()
	sFBM_MSDs[:, j] = [k[0] for k in output]

ls = np.zeros([len(penalties)])

for i in range(ls.size):
	ls[i] = np.sum(np.square(sFBM_MSDs[:, i] - MD_MSDs))

print('Optimal penalty: %s' % penalties[np.argmin(ls)])
