#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sqlite3 as sql

path = "/home/bcoscia/PycharmProjects/LLC_Membranes/LLC_Membranes/analysis/"
connection = sql.connect('%s/size.db' % path)
crsr = connection.cursor()
wt = 10 

simple_alcohols = True 
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

command = "SELECT name, MD_TAMSD, MD_TAMSD_CI_lower, MD_TAMSD_CI_upper, sim_length FROM msd WHERE wt_water = %d and name != 'HOH'" % wt

if wt == 10:
	command += "and penalty = 0.25"

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

# Theoretical Calcuations
water_size = 0.151391485766 # from ACH cubic box. std: 0.000407801323403

r_theoretical = np.linspace(0.25, 0.68, 100)
#r_theoretical = np.linspace(min(r), 10, 1000)

# Corrected version first
alpha = water_size / r_theoretical # radius of solute (residue) to solvent (water) 
f = ((3 * alpha / 2) + (1 / (1 + alpha)))**-1  # f as function of alpha
alpha_max = water_size / r[r_order][np.argmax(msd)]  # alpha for solute with max msd
fmax = ((3 * alpha_max / 2) + (1 / (1 + alpha_max)))**-1  # f for solute with max msd
a = fmax * np.max(msd) * r[r_order][np.argmax(msd)]  # make the curve pass through the max msd point

# Make Stokes-Einstein intersect with correction at r_intersect
r_intersect = 5 # where we want stokes-einstein and corrected stoke-einstein to intersect

alpha_intersect = water_size / r_intersect
f_intersect = ((3 * alpha_intersect / 2) + (1 / (1 + alpha_intersect)))**-1
y_intersect = a / (f_intersect * r_intersect)
y_stokes = y_intersect * r_intersect / r_theoretical
y_gierer_wirtz = a / (r_theoretical * f)

plt.plot(r_theoretical, y_stokes, '--', color='black', label='Stokes-Einstein')
plt.plot(r_theoretical, y_gierer_wirtz, '--', color='blue', label='Gierer and Wirtz correction')

# Plot formatting
plt.xlabel('Radius (nm)', fontsize=14)
plt.ylabel('MSD ($nm^2$)', fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.legend(fontsize=14)
if restrict:
	plt.savefig('msd_radius_simple_alcohols_%dwt.pdf' % wt)
else:
	plt.savefig('msd_radius_%dwt.pdf' % wt)
plt.show()
