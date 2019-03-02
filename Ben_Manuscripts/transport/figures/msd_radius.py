#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql

path = "/home/bcoscia/PycharmProjects/LLC_Membranes/LLC_Membranes/analysis/"
connection = sql.connect('%s/size.db' % path)
crsr = connection.cursor()

simple_alcohols = False 
if simple_alcohols:
	restrict = ['MET', 'ETH', 'PR', 'BUT']
else:
	restrict = []

command = "SELECT name, end_to_end, end_to_end_std FROM size WHERE nsolute = 6"

if restrict:
	command += ' and ('
	for i, n in enumerate(restrict):
		command += " name = '%s'" % n 
		if i < len(restrict) - 1:
			command += " or"
	command += ')'

radii = crsr.execute(command).fetchall()

rnames = [i[0] for i in radii]
r = np.array([i[1] for i in radii])
r_std = np.array([i[2] for i in radii])

connection = sql.connect('%s/../timeseries/msd.db' % path)
crsr = connection.cursor()

command = "SELECT name, MD_TAMSD, MD_TAMSD_CI_lower, MD_TAMSD_CI_upper, sim_length FROM msd WHERE penalty = 0.25"

if restrict:
	command += ' and ('
	for i, n in enumerate(restrict):
		command += " name = '%s'" % n 
		if i < len(restrict) - 1:
			command += " or"
	command += ')'

msds = crsr.execute(command).fetchall()

names = np.array([i[0] for i in msds])
msd = np.array([i[1] for i in msds])
msd_lower = -np.array([i[2] for i in msds]) + msd
msd_upper = np.array([i[3] for i in msds]) - msd
sim_length = np.array([i[4] for i in msds])

msd *= (400 / (sim_length * 0.4))

ordered = np.argsort(msd)

names = names[ordered].tolist()
msd = msd[ordered]
msd_lower = msd_lower[ordered]
msd_upper = msd_upper[ordered]
sim_length = sim_length[ordered]

r_order = np.array([rnames.index(i) for i in names])

#plt.scatter(msd, r[r_order])
plt.errorbar(r[r_order], msd, yerr=(msd_lower, msd_upper), xerr=r_std[r_order], fmt='o', elinewidth=0.5, ecolor='black')
plt.xlabel('Radius (nm)', fontsize=14)
plt.ylabel('MSD ($nm^2$)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
if restrict:
	plt.savefig('msd_radius_simple_alcohols.pdf')
else:
	plt.savefig('msd_radius.pdf')
plt.show()
